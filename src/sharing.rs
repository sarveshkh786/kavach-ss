//! Kavach-SS: Cauchy-MDS threshold sharing with an information-theoretic
//! MAC that lets the combiner *validate* the reconstructed secret.
//!
//! For a secret split into L limbs s_1..s_L (56 bits each) the dealer samples
//! a MAC key alpha and shares the 2L+1 field values
//!     (alpha, s_1, alpha*s_1, ..., s_L, alpha*s_L).
//! Each value w is shared with an independent random vector v in F^t chosen
//! uniformly subject to C_0 . v = w; share i receives C_i . v.

use crate::cauchy;
use crate::field::{zeroize, zeroize_bytes, Fp, P};
use std::fmt;
use std::fs::File;
use std::io::Read;

pub const LIMB_BYTES: usize = 7; // 56-bit limbs < p
pub const MAGIC: &str = "kavach1";
pub const MAX_N: usize = 255;
pub const MAX_SECRET: usize = 1 << 16;

#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    BadParams(&'static str),
    NotEnoughShares {
        have: usize,
        need: usize,
    },
    Mismatched,
    DuplicateIndex(u16),
    /// Redundant share disagrees with the others (detected via MDS property).
    Inconsistent(u16),
    /// Information-theoretic MAC check failed: tampering detected.
    MacFailure,
    Parse(&'static str),
    Rng,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::BadParams(s) => write!(f, "bad parameters: {s}"),
            Error::NotEnoughShares { have, need } => {
                write!(f, "not enough shares: have {have}, need {need}")
            }
            Error::Mismatched => write!(f, "shares come from different dealings"),
            Error::DuplicateIndex(i) => write!(f, "duplicate share index {i}"),
            Error::Inconsistent(i) => write!(
                f,
                "share set is INCONSISTENT (mismatch detected at share #{i}): a share was tampered or corrupted"
            ),
            Error::MacFailure => write!(
                f,
                "MAC validation FAILED: reconstructed secret is not authentic"
            ),
            Error::Parse(s) => write!(f, "parse error: {s}"),
            Error::Rng => write!(f, "OS randomness unavailable"),
        }
    }
}

/// OS CSPRNG (/dev/urandom); no third-party code.
pub struct Rng {
    f: File,
}

impl Rng {
    pub fn new() -> Result<Rng, Error> {
        Ok(Rng {
            f: File::open("/dev/urandom").map_err(|_| Error::Rng)?,
        })
    }
    pub fn bytes(&mut self, buf: &mut [u8]) -> Result<(), Error> {
        self.f.read_exact(buf).map_err(|_| Error::Rng)
    }
    /// Uniform field element by rejection sampling on 61-bit draws.
    /// (Rejection depends only on fresh randomness, never on the secret.)
    pub fn fp(&mut self) -> Result<Fp, Error> {
        loop {
            let mut b = [0u8; 8];
            self.bytes(&mut b)?;
            let v = u64::from_le_bytes(b) & P;
            zeroize_bytes(&mut b);
            if v != P {
                return Ok(Fp(v));
            }
        }
    }
    pub fn fp_nonzero(&mut self) -> Result<Fp, Error> {
        loop {
            let v = self.fp()?;
            if v.ct_eq(Fp::ZERO) == 0 {
                return Ok(v);
            }
        }
    }
}

/// One share. Index, t, n, length and group id are public metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Share {
    pub index: u16,
    pub t: u16,
    pub n: u16,
    pub secret_len: u32,
    pub group: [u8; 16],
    pub values: Vec<Fp>,
}

impl Drop for Share {
    fn drop(&mut self) {
        zeroize(&mut self.values);
    }
}

/// A byte buffer that is wiped when dropped.
pub struct SecretBytes(pub Vec<u8>);
impl Drop for SecretBytes {
    fn drop(&mut self) {
        zeroize_bytes(&mut self.0);
    }
}

fn limbs_of(secret: &[u8]) -> Vec<Fp> {
    secret
        .chunks(LIMB_BYTES)
        .map(|c| {
            let mut b = [0u8; 8];
            b[..c.len()].copy_from_slice(c);
            let v = Fp(u64::from_le_bytes(b));
            zeroize_bytes(&mut b);
            v
        })
        .collect()
}

/// Share a single field value w. Returns the n share-components.
fn share_value(
    w: Fp,
    t: usize,
    n: usize,
    c0: &[Fp],
    c00_inv: Fp,
    rng: &mut Rng,
) -> Result<Vec<Fp>, Error> {
    let mut v = vec![Fp::ZERO; t];
    let mut acc = Fp::ZERO;
    for j in 1..t {
        v[j] = rng.fp()?;
        acc = acc.add(c0[j].mul(v[j]));
    }
    v[0] = w.sub(acc).mul(c00_inv);
    let out = (1..=n)
        .map(|i| {
            let r = cauchy::row(i, t);
            r.iter()
                .zip(&v)
                .fold(Fp::ZERO, |a, (c, x)| a.add(c.mul(*x)))
        })
        .collect();
    zeroize(&mut v);
    Ok(out)
}

/// Split `secret` into n shares, any t of which reconstruct it.
pub fn split(secret: &[u8], t: usize, n: usize, rng: &mut Rng) -> Result<Vec<Share>, Error> {
    if t == 0 || t > n {
        return Err(Error::BadParams("need 1 <= t <= n"));
    }
    if n > MAX_N {
        return Err(Error::BadParams("n must be <= 255"));
    }
    if secret.is_empty() || secret.len() > MAX_SECRET {
        return Err(Error::BadParams("secret must be 1..=65536 bytes"));
    }
    let c0 = cauchy::row(0, t);
    let c00_inv = c0[0].inv();
    let mut group = [0u8; 16];
    rng.bytes(&mut group)?;

    let mut limbs = limbs_of(secret);
    let mut alpha = rng.fp_nonzero()?;
    let mut plain: Vec<Fp> = Vec::with_capacity(2 * limbs.len() + 1);
    plain.push(alpha);
    for s in &limbs {
        plain.push(*s);
        plain.push(alpha.mul(*s));
    }

    let mut shares: Vec<Share> = (1..=n)
        .map(|i| Share {
            index: i as u16,
            t: t as u16,
            n: n as u16,
            secret_len: secret.len() as u32,
            group,
            values: Vec::with_capacity(plain.len()),
        })
        .collect();
    for w in &plain {
        let comps = share_value(*w, t, n, &c0, c00_inv, rng)?;
        for (sh, c) in shares.iter_mut().zip(comps) {
            sh.values.push(c);
        }
    }
    zeroize(&mut plain);
    zeroize(&mut limbs);
    alpha = Fp::ZERO;
    let _ = alpha;
    Ok(shares)
}

fn check_set(shares: &[Share]) -> Result<(usize, usize), Error> {
    let first = shares
        .first()
        .ok_or(Error::NotEnoughShares { have: 0, need: 1 })?;
    let (t, n) = (first.t as usize, first.n as usize);
    let width = 2 * first.secret_len.div_ceil(LIMB_BYTES as u32) as usize + 1;
    for s in shares {
        if s.t != first.t
            || s.n != first.n
            || s.secret_len != first.secret_len
            || s.group != first.group
        {
            return Err(Error::Mismatched);
        }
        if s.values.len() != width || s.index == 0 || s.index as usize > n {
            return Err(Error::Parse("malformed share"));
        }
    }
    for (a, s) in shares.iter().enumerate() {
        if shares[..a].iter().any(|o| o.index == s.index) {
            return Err(Error::DuplicateIndex(s.index));
        }
    }
    if shares.len() < t {
        return Err(Error::NotEnoughShares {
            have: shares.len(),
            need: t,
        });
    }
    Ok((t, n))
}

/// Reconstruct and validate. Uses the first t shares; any extra shares are
/// checked for consistency (error detection from the MDS code), then the MAC
/// tag alpha*s_l is verified for every limb in constant time.
pub fn combine(shares: &[Share]) -> Result<SecretBytes, Error> {
    let (t, _n) = check_set(shares)?;
    let base: Vec<usize> = shares[..t].iter().map(|s| s.index as usize).collect();
    let width = shares[0].values.len();

    // Secret recombination vector lambda = C_0 . C_I^{-1}
    let lambda = cauchy::recombination(&base, 0, t).ok_or(Error::BadParams("singular"))?;

    // Consistency of redundant shares: share_e must equal mu_e . shares_I
    for extra in &shares[t..] {
        let mu = cauchy::recombination(&base, extra.index as usize, t)
            .ok_or(Error::BadParams("singular"))?;
        let mut ok = 1u64;
        for w in 0..width {
            let pred = mu
                .iter()
                .zip(&shares[..t])
                .fold(Fp::ZERO, |a, (m, s)| a.add(m.mul(s.values[w])));
            ok &= pred.ct_eq(extra.values[w]);
        }
        if ok != 1 {
            return Err(Error::Inconsistent(extra.index));
        }
    }

    let mut vals: Vec<Fp> = (0..width)
        .map(|w| {
            lambda
                .iter()
                .zip(&shares[..t])
                .fold(Fp::ZERO, |a, (l, s)| a.add(l.mul(s.values[w])))
        })
        .collect();

    // Constant-time MAC verification over all limbs, plus range check
    let alpha = vals[0];
    let mut ok = 1u64;
    let limbs = (width - 1) / 2;
    for l in 0..limbs {
        let s = vals[1 + 2 * l];
        let tag = vals[2 + 2 * l];
        ok &= alpha.mul(s).ct_eq(tag);
        ok &= 1 ^ ((s.value() >> 56) != 0) as u64; // limb must fit 56 bits
    }
    ok &= 1 ^ alpha.ct_eq(Fp::ZERO);
    if ok != 1 {
        zeroize(&mut vals);
        return Err(Error::MacFailure);
    }

    let len = shares[0].secret_len as usize;
    let mut out = SecretBytes(Vec::with_capacity(limbs * LIMB_BYTES));
    for l in 0..limbs {
        out.0
            .extend_from_slice(&vals[1 + 2 * l].value().to_le_bytes()[..LIMB_BYTES]);
    }
    zeroize(&mut vals);
    // Trailing padding must be zero
    let pad_ok = out.0[len..].iter().fold(0u8, |a, b| a | b) == 0;
    if !pad_ok {
        return Err(Error::MacFailure);
    }
    out.0.truncate(len);
    Ok(out)
}

/// Proactive refresh: re-randomise all shares by adding a fresh sharing of
/// zero to every component. The secret and MAC relation are unchanged, but
/// shares stolen before the refresh become useless when combined with shares
/// taken after it. Requires the full share set in this simulation; in a
/// deployment each holder deals a zero-sharing and everyone sums.
pub fn refresh(shares: &mut [Share], rng: &mut Rng) -> Result<(), Error> {
    let (t, n) = check_set(shares)?;
    if shares.len() != n {
        return Err(Error::BadParams("refresh needs all n shares"));
    }
    let c0 = cauchy::row(0, t);
    let c00_inv = c0[0].inv();
    let width = shares[0].values.len();
    let mut newgroup = [0u8; 16];
    rng.bytes(&mut newgroup)?;
    for w in 0..width {
        let z = share_value(Fp::ZERO, t, n, &c0, c00_inv, rng)?;
        for sh in shares.iter_mut() {
            let k = sh.index as usize - 1;
            sh.values[w] = sh.values[w].add(z[k]);
        }
    }
    for sh in shares.iter_mut() {
        sh.group = newgroup; // epoch change: old and new shares refuse to mix
    }
    Ok(())
}

// ---------- serialisation (text, line-safe) ----------

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}
fn unhex(s: &str) -> Result<Vec<u8>, Error> {
    if s.len() % 2 != 0 {
        return Err(Error::Parse("odd hex"));
    }
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).map_err(|_| Error::Parse("bad hex")))
        .collect()
}

impl Share {
    pub fn encode(&self) -> String {
        let mut b = Vec::new();
        b.extend_from_slice(&self.index.to_le_bytes());
        b.extend_from_slice(&self.t.to_le_bytes());
        b.extend_from_slice(&self.n.to_le_bytes());
        b.extend_from_slice(&self.secret_len.to_le_bytes());
        b.extend_from_slice(&self.group);
        for v in &self.values {
            b.extend_from_slice(&v.value().to_le_bytes());
        }
        let s = format!("{MAGIC}-{}of{}-{}", self.t, self.n, hex(&b));
        zeroize_bytes(&mut b);
        s
    }

    pub fn decode(s: &str) -> Result<Share, Error> {
        let s = s.trim();
        let mut parts = s.splitn(3, '-');
        if parts.next() != Some(MAGIC) {
            return Err(Error::Parse("not a kavach share"));
        }
        let _hdr = parts.next().ok_or(Error::Parse("missing header"))?;
        let mut b = unhex(parts.next().ok_or(Error::Parse("missing body"))?)?;
        if b.len() < 26 || (b.len() - 26) % 8 != 0 {
            return Err(Error::Parse("bad length"));
        }
        let u16at = |o: usize| u16::from_le_bytes([b[o], b[o + 1]]);
        let mut group = [0u8; 16];
        group.copy_from_slice(&b[10..26]);
        let mut values = Vec::new();
        for c in b[26..].chunks(8) {
            let v = u64::from_le_bytes(c.try_into().unwrap());
            if v >= P {
                return Err(Error::Parse("non-canonical field element"));
            }
            values.push(Fp(v));
        }
        let sh = Share {
            index: u16at(0),
            t: u16at(2),
            n: u16at(4),
            secret_len: u32::from_le_bytes([b[6], b[7], b[8], b[9]]),
            group,
            values,
        };
        zeroize_bytes(&mut b);
        Ok(sh)
    }
}
