# Kavach-SS — validated threshold sharing for master keys

**Kavach** (Hindi: *armour*) splits a master key into `n` shares so that:

* any **t** shares reconstruct it,
* any **t − 1** shares reveal **nothing**: perfect, information-theoretic secrecy that holds against quantum or unbounded adversaries,
* the reconstructor can **validate** the result. A single field-element MAC key shared alongside the secret catches any tampering by up to t − 1 liars with probability ≥ 1 − 2⁻⁶⁰. With redundant shares, detection is certain and the bad share is named,
* shares can be **proactively refreshed**, so shares stolen in an earlier epoch become useless.

The maths is number theory and linear algebra only: the Mersenne prime field **GF(2⁶¹ − 1)** and an **extended Cauchy matrix**, whose every square minor is non-zero. Reconstruction and privacy both follow from that one determinant identity. No existing cipher touches the plaintext.

📄 **Full write-up with proofs:** [`paper/kavach.pdf`](paper/kavach.pdf) (LaTeX source: [`paper/kavach.tex`](paper/kavach.tex))

## Build & run

```bash
cargo build --release          # zero dependencies
cargo test --release           # 11 security tests
./target/release/kavach demo   # end-to-end demo with attack scenarios
```

```
kavach genkey  -o master.key
kavach split   -t 3 -n 5 -i master.key -o shares/
kavach combine -o rec.key shares/share-2.kv shares/share-4.kv shares/share-5.kv
kavach verify  shares/*.kv          # validate without writing the key
kavach refresh shares/*.kv          # proactive re-randomisation (new epoch)
kavach audit   -t 3 -n 5            # check every reconstruction + privacy minor
kavach ctcheck                      # dudect-style timing-leak test
./scripts/asm_audit.sh              # count branches in compiled field arithmetic
```

## Demo output

```
[1] combine shares {1,3,5}            -> recovered = original: true
[2] combine {4,2,1} + redundant #3     -> consistent & recovered: true
[3] attacker holds only 2 shares         -> not enough shares: have 2, need 3
[4] malicious holder flips share #2 bits -> MAC validation FAILED: reconstructed secret is not authentic
[5] same attack, 4 shares supplied       -> share #4 is inconsistent with the others (tampered or corrupted)
[6] proactive refresh, then combine      -> still recovers key: true
[7] share stolen BEFORE refresh + 2 new  -> shares come from different dealings
```

## Side-channel hardening

| Primitive (rustc -O, x86-64) | instructions | conditional branches |
|---|---|---|
| `Fp::mul` | 22 | 0 |
| `Fp::add` / `Fp::sub` | 10 | 0 |
| `Fp::ct_eq` (MAC compare) | 4 | 0 |
| Mersenne reduction | 12 | 0 |
| `Fp::pow` | 46 | 1 (fixed 64-iteration loop counter) |

* Branch-free conditional subtraction; fixed-schedule exponentiation; the MAC is checked across all limbs before one accept/reject decision.
* Only public data (share indices, lengths) steers control flow.
* All secret buffers are zeroised with `write_volatile` and a compiler fence; output files are created with mode 0600.
* No third-party crates; randomness comes from the kernel CSPRNG.

## Layout

```
src/field.rs     constant-time GF(2^61-1)
src/cauchy.rs    extended Cauchy generator, public linear algebra
src/sharing.rs   deal / combine+validate / refresh / encoding
src/ctcheck.rs   dudect-style Welch t-test
src/main.rs      CLI
tests/           security test-suite
paper/           LaTeX write-up
```

License: MIT
