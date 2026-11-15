use kavach::field::Fp;
use kavach::sharing::*;

fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn go(n: usize, k: usize, s: usize, c: &mut Vec<usize>, o: &mut Vec<Vec<usize>>) {
        if c.len() == k {
            o.push(c.clone());
            return;
        }
        for i in s..n {
            c.push(i);
            go(n, k, i + 1, c, o);
            c.pop();
        }
    }
    let mut o = vec![];
    go(n, k, 0, &mut vec![], &mut o);
    o
}

#[test]
fn every_t_subset_reconstructs() {
    let mut rng = Rng::new().unwrap();
    for (t, n) in [(1, 1), (1, 4), (2, 3), (3, 5), (4, 7), (5, 5), (6, 9)] {
        for len in [1usize, 7, 31, 32, 64] {
            let mut key = vec![0u8; len];
            rng.bytes(&mut key).unwrap();
            let sh = split(&key, t, n, &mut rng).unwrap();
            for s in subsets(n, t) {
                let pick: Vec<Share> = s.iter().map(|&i| sh[i].clone()).collect();
                assert_eq!(combine(&pick).unwrap().0, key, "t={t} n={n} {s:?}");
            }
            // all n together (redundancy check path)
            assert_eq!(combine(&sh).unwrap().0, key);
        }
    }
}

#[test]
fn below_threshold_is_rejected() {
    let mut rng = Rng::new().unwrap();
    let sh = split(b"master-key-0123456789abcdef!!!!!", 4, 6, &mut rng).unwrap();
    assert_eq!(
        combine(&sh[..3]).err(),
        Some(Error::NotEnoughShares { have: 3, need: 4 })
    );
}

/// Any single-component tamper by one holder is caught (by MAC or redundancy).
#[test]
fn tampering_is_detected_every_position() {
    let mut rng = Rng::new().unwrap();
    let key = [0x42u8; 32];
    let sh = split(&key, 3, 5, &mut rng).unwrap();
    for who in 0..3 {
        for pos in 0..sh[0].values.len() {
            for delta in [1u64, 2, 1 << 40, (1 << 61) - 2] {
                let mut s = vec![sh[0].clone(), sh[1].clone(), sh[2].clone()];
                s[who].values[pos] = s[who].values[pos].add(Fp::new(delta));
                assert_eq!(
                    combine(&s).err(),
                    Some(Error::MacFailure),
                    "who={who} pos={pos}"
                );
                s.push(sh[3].clone());
                assert!(matches!(combine(&s).err(), Some(Error::Inconsistent(_))));
            }
        }
    }
}

/// Adaptive attacker: holder of share #1 tries to shift the secret by a chosen
/// offset AND fix up the tag. Without alpha the success rate must be ~1/p.
#[test]
fn targeted_forgery_fails() {
    let mut rng = Rng::new().unwrap();
    let mut wins = 0;
    for _ in 0..500 {
        let sh = split(&[7u8; 14], 2, 3, &mut rng).unwrap();
        let mut s = vec![sh[0].clone(), sh[1].clone()];
        // shift limb 0 of the secret by +1 at reconstruction: adjust component 1,
        // and guess a tag correction (attacker does not know alpha)
        let lam = kavach::cauchy::recombination(&[1, 2], 0, 2).unwrap();
        let fix = lam[0].inv();
        s[0].values[1] = s[0].values[1].add(fix);
        s[0].values[2] = s[0].values[2].add(rng.fp().unwrap().mul(fix));
        if combine(&s).is_ok() {
            wins += 1;
        }
    }
    assert_eq!(wins, 0);
}

#[test]
fn mixing_dealings_is_rejected() {
    let mut rng = Rng::new().unwrap();
    let a = split(&[1u8; 32], 2, 3, &mut rng).unwrap();
    let b = split(&[1u8; 32], 2, 3, &mut rng).unwrap();
    assert_eq!(
        combine(&[a[0].clone(), b[1].clone()]).err(),
        Some(Error::Mismatched)
    );
}

#[test]
fn refresh_preserves_secret_and_kills_old_shares() {
    let mut rng = Rng::new().unwrap();
    let key = b"0123456789abcdef0123456789abcdef";
    let mut sh = split(key, 3, 5, &mut rng).unwrap();
    let old = sh.clone();
    refresh(&mut sh, &mut rng).unwrap();
    assert_eq!(combine(&sh[2..5]).unwrap().0, key);
    for i in 0..5 {
        assert_ne!(old[i].values, sh[i].values);
    }
    // even if an attacker forges the group id, old+new shares fail the MAC
    let mut stale = old[0].clone();
    stale.group = sh[1].group;
    assert!(combine(&[stale, sh[1].clone(), sh[2].clone()]).is_err());
}

#[test]
fn encode_decode_roundtrip_and_reject_garbage() {
    let mut rng = Rng::new().unwrap();
    let sh = split(&[9u8; 32], 2, 2, &mut rng).unwrap();
    for s in &sh {
        assert_eq!(&Share::decode(&s.encode()).unwrap(), s);
    }
    assert!(Share::decode("hello").is_err());
    assert!(Share::decode("kavach1-2of2-zz").is_err());
    let mut e = sh[0].encode();
    e.push_str("ffffffffffffffff"); // non-canonical element
    assert!(Share::decode(&e).is_err());
}

/// Empirical privacy check: with t-1 shares fixed, the distribution of a
/// share component is independent of the secret (chi-square on low bits).
#[test]
fn t_minus_1_shares_look_uniform() {
    let mut rng = Rng::new().unwrap();
    let mut counts = [[0u32; 16]; 2];
    for (c, key) in [[0u8; 7], [0xFFu8; 7]].iter().enumerate() {
        for _ in 0..4000 {
            let sh = split(key, 3, 5, &mut rng).unwrap();
            counts[c][(sh[0].values[1].value() & 15) as usize] += 1;
        }
    }
    for c in 0..2 {
        let chi: f64 = counts[c]
            .iter()
            .map(|&o| {
                let e = 250.0;
                (o as f64 - e).powi(2) / e
            })
            .sum();
        assert!(chi < 37.7, "chi2 = {chi}"); // p = 0.001, 15 dof
    }
}
