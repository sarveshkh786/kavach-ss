use kavach::sharing::{combine, refresh, split, Rng, SecretBytes, Share};
use kavach::{cauchy, ctcheck};
use std::fs;
use std::io::Write;
use std::process::exit;

const USAGE: &str = "\
Kavach-SS — validated threshold secret sharing for master keys (Rust, zero dependencies)

USAGE:
  kavach split   -t <T> -n <N> -i <secret-file> -o <out-dir>
  kavach combine -o <out-file> <share-file>...
  kavach verify  <share-file>...          check shares + MAC without writing the secret
  kavach refresh <share-file>...          proactively re-randomise ALL n shares in place
  kavach genkey  -o <file> [-b bytes]     create a random master key (default 32 bytes)
  kavach audit   -t <T> -n <N>            prove MDS + privacy conditions for these params
  kavach ctcheck [-s samples]             dudect timing-leak test of field arithmetic
  kavach demo                             end-to-end demonstration with attack scenarios
";

fn die(msg: impl std::fmt::Display) -> ! {
    eprintln!("error: {msg}");
    exit(1)
}

fn flag(args: &[String], f: &str) -> Option<String> {
    args.iter()
        .position(|a| a == f)
        .and_then(|i| args.get(i + 1).cloned())
}

fn positional(args: &[String]) -> Vec<String> {
    let mut out = vec![];
    let mut i = 0;
    while i < args.len() {
        if args[i].starts_with('-') {
            i += 2;
        } else {
            out.push(args[i].clone());
            i += 1;
        }
    }
    out
}

fn load_shares(paths: &[String]) -> Vec<Share> {
    paths
        .iter()
        .map(|p| {
            let s = fs::read_to_string(p).unwrap_or_else(|e| die(format!("{p}: {e}")));
            Share::decode(&s).unwrap_or_else(|e| die(format!("{p}: {e}")))
        })
        .collect()
}

fn write_private(path: &str, data: &[u8]) {
    let mut opts = fs::OpenOptions::new();
    opts.write(true).create(true).truncate(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        opts.mode(0o600);
    }
    let mut f = opts
        .open(path)
        .unwrap_or_else(|e| die(format!("{path}: {e}")));
    f.write_all(data).unwrap_or_else(|e| die(e));
}

fn audit(t: usize, n: usize) -> bool {
    fn subsets(
        n: usize,
        k: usize,
        start: usize,
        cur: &mut Vec<usize>,
        f: &mut dyn FnMut(&[usize]) -> bool,
    ) -> bool {
        if cur.len() == k {
            return f(cur);
        }
        for i in start..=n {
            cur.push(i);
            if !subsets(n, k, i + 1, cur, f) {
                return false;
            }
            cur.pop();
        }
        true
    }
    let mut count = 0u64;
    let recon = subsets(n, t, 1, &mut vec![], &mut |s| {
        count += 1;
        cauchy::is_nonsingular(&s.iter().map(|&i| cauchy::row(i, t)).collect::<Vec<_>>())
    });
    println!("  reconstruction: {count} subsets of size t={t} -> all invertible: {recon}");
    let mut count2 = 0u64;
    let privacy = subsets(n, t - 1, 1, &mut vec![], &mut |s| {
        count2 += 1;
        let mut rows = vec![cauchy::row(0, t)];
        rows.extend(s.iter().map(|&i| cauchy::row(i, t)));
        cauchy::is_nonsingular(&rows)
    });
    println!(
        "  perfect privacy: {count2} coalitions of size t-1={} -> [C_0; C_J] invertible: {privacy}",
        t - 1
    );
    recon && privacy
}

fn demo() {
    let mut rng = Rng::new().unwrap_or_else(|e| die(e));
    println!("== Kavach-SS demo ==\n");
    let mut key = SecretBytes(vec![0u8; 32]);
    rng.bytes(&mut key.0).unwrap();
    let fp: String = key.0[..4].iter().map(|b| format!("{b:02x}")).collect();
    println!("[dealer] generated 256-bit master key (fingerprint {fp}...)");
    let mut shares = split(&key.0, 3, 5, &mut rng).unwrap();
    println!("[dealer] split into n=5 shares, threshold t=3; each share = {} field elements over GF(2^61-1)", shares[0].values.len());
    for s in &shares {
        let e = s.encode();
        println!(
            "   share #{}: {}...{}",
            s.index,
            &e[..34],
            &e[e.len() - 8..]
        );
    }

    let r = combine(&[shares[0].clone(), shares[2].clone(), shares[4].clone()]).unwrap();
    println!(
        "\n[1] combine shares {{1,3,5}}            -> recovered = original: {}",
        r.0 == key.0
    );
    let r = combine(&[
        shares[3].clone(),
        shares[1].clone(),
        shares[0].clone(),
        shares[2].clone(),
    ])
    .unwrap();
    println!(
        "[2] combine {{4,2,1}} + redundant #3     -> consistent & recovered: {}",
        r.0 == key.0
    );

    let e = combine(&shares[..2].to_vec()).err().unwrap();
    println!("[3] attacker holds only 2 shares         -> {e}");

    let mut bad = shares[1].clone();
    bad.values[3] = bad.values[3].add(kavach::field::Fp::new(1));
    let e = combine(&[shares[0].clone(), bad.clone(), shares[2].clone()])
        .err()
        .unwrap();
    println!("[4] malicious holder flips share #2 bits -> {e}");
    let e = combine(&[shares[0].clone(), bad, shares[2].clone(), shares[3].clone()])
        .err()
        .unwrap();
    println!("[5] same attack, 4 shares supplied       -> {e}");

    let stolen = shares[0].clone();
    refresh(&mut shares, &mut rng).unwrap();
    let r = combine(&[shares[1].clone(), shares[3].clone(), shares[4].clone()]).unwrap();
    println!(
        "[6] proactive refresh, then combine      -> still recovers key: {}",
        r.0 == key.0
    );
    let e = combine(&[stolen, shares[1].clone(), shares[2].clone()])
        .err()
        .unwrap();
    println!("[7] share stolen BEFORE refresh + 2 new  -> {e}");
    println!("\nAll secret buffers are zeroised on drop.");
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = argv.first() else {
        print!("{USAGE}");
        exit(2)
    };
    let args = &argv[1..];
    match cmd.as_str() {
        "split" => {
            let t: usize = flag(args, "-t")
                .and_then(|v| v.parse().ok())
                .unwrap_or_else(|| die("missing -t"));
            let n: usize = flag(args, "-n")
                .and_then(|v| v.parse().ok())
                .unwrap_or_else(|| die("missing -n"));
            let inp = flag(args, "-i").unwrap_or_else(|| die("missing -i"));
            let out = flag(args, "-o").unwrap_or_else(|| die("missing -o"));
            let secret = SecretBytes(fs::read(&inp).unwrap_or_else(|e| die(format!("{inp}: {e}"))));
            let mut rng = Rng::new().unwrap_or_else(|e| die(e));
            let shares = split(&secret.0, t, n, &mut rng).unwrap_or_else(|e| die(e));
            fs::create_dir_all(&out).unwrap_or_else(|e| die(e));
            for s in &shares {
                let p = format!("{out}/share-{}.kv", s.index);
                write_private(&p, s.encode().as_bytes());
                println!("wrote {p}");
            }
            println!(
                "{n} shares written; any {t} reconstruct the secret, any {} reveal nothing.",
                t - 1
            );
        }
        "combine" => {
            let out = flag(args, "-o").unwrap_or_else(|| die("missing -o"));
            let shares = load_shares(&positional(args));
            let secret = combine(&shares).unwrap_or_else(|e| die(e));
            write_private(&out, &secret.0);
            println!(
                "MAC valid. Secret ({} bytes) written to {out} (mode 0600).",
                secret.0.len()
            );
        }
        "verify" => {
            let shares = load_shares(&positional(args));
            match combine(&shares) {
                Ok(s) => println!(
                    "OK: {} shares consistent, MAC valid, secret length {} bytes.",
                    shares.len(),
                    s.0.len()
                ),
                Err(e) => die(e),
            }
        }
        "refresh" => {
            let paths = positional(args);
            let mut shares = load_shares(&paths);
            let mut rng = Rng::new().unwrap_or_else(|e| die(e));
            refresh(&mut shares, &mut rng).unwrap_or_else(|e| die(e));
            for (p, s) in paths.iter().zip(&shares) {
                write_private(p, s.encode().as_bytes());
            }
            println!(
                "Refreshed {} shares (new epoch). Old shares can no longer be mixed in.",
                shares.len()
            );
        }
        "genkey" => {
            let out = flag(args, "-o").unwrap_or_else(|| die("missing -o"));
            let b: usize = flag(args, "-b").and_then(|v| v.parse().ok()).unwrap_or(32);
            let mut k = SecretBytes(vec![0u8; b]);
            Rng::new()
                .unwrap_or_else(|e| die(e))
                .bytes(&mut k.0)
                .unwrap_or_else(|e| die(e));
            write_private(&out, &k.0);
            println!("wrote {b}-byte random key to {out}");
        }
        "audit" => {
            let t: usize = flag(args, "-t")
                .and_then(|v| v.parse().ok())
                .unwrap_or_else(|| die("missing -t"));
            let n: usize = flag(args, "-n")
                .and_then(|v| v.parse().ok())
                .unwrap_or_else(|| die("missing -n"));
            if t == 0 || t > n || n > 30 {
                die("need 1 <= t <= n <= 30 for exhaustive audit");
            }
            println!("Auditing Cauchy generator for t={t}, n={n}:");
            if audit(t, n) {
                println!("PASS");
            } else {
                die("FAIL");
            }
        }
        "ctcheck" => {
            let s: usize = flag(args, "-s")
                .and_then(|v| v.parse().ok())
                .unwrap_or(200_000);
            println!("dudect-style Welch t-test, fixed(0) vs random secret class; |t| < 4.5 => no leak detected");
            let mut all_ok = true;
            for r in ctcheck::run_all(s) {
                let ok = r.t.abs() < 4.5;
                all_ok &= ok;
                println!(
                    "  {:<30} n={:<8} t = {:>7.3}  {}",
                    r.name,
                    r.samples,
                    r.t,
                    if ok { "PASS" } else { "LEAK?" }
                );
            }
            if !all_ok {
                exit(1)
            }
        }
        "demo" => demo(),
        _ => {
            print!("{USAGE}");
            exit(2)
        }
    }
}
