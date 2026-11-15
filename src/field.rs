//! Constant-time arithmetic in the Mersenne prime field GF(p), p = 2^61 - 1.
//!
//! Every operation on secret-bearing elements is branch-free and uses no
//! secret-dependent memory indexing. Reduction exploits 2^61 = 1 (mod p).

use core::sync::atomic::{compiler_fence, Ordering};

/// The Mersenne prime 2^61 - 1.
pub const P: u64 = (1u64 << 61) - 1;

/// A field element. Invariant: `0 <= self.0 < P`.
#[derive(Clone, Copy, Debug, Default)]
pub struct Fp(pub(crate) u64);

/// Branch-free conditional subtraction of P: returns r - P if r >= P, else r.
/// Valid for r < 2^63.
#[inline(always)]
fn ct_sub_p(r: u64) -> u64 {
    let s = r.wrapping_sub(P);
    // If r < P the subtraction wraps and the top bit of s is 1.
    let keep_r = 0u64.wrapping_sub(s >> 63); // all-ones iff r < P
    (s & !keep_r) | (r & keep_r)
}

impl Fp {
    pub const ZERO: Fp = Fp(0);
    pub const ONE: Fp = Fp(1);

    /// Reduce any u64 into the field (constant time).
    #[inline(always)]
    pub fn new(x: u64) -> Fp {
        let r = (x & P) + (x >> 61); // < 2^61 + 8
        Fp(ct_sub_p(r))
    }

    #[inline(always)]
    pub fn value(self) -> u64 {
        self.0
    }

    #[inline(always)]
    pub fn add(self, o: Fp) -> Fp {
        Fp(ct_sub_p(self.0 + o.0))
    }

    #[inline(always)]
    pub fn sub(self, o: Fp) -> Fp {
        Fp(ct_sub_p(self.0 + P - o.0))
    }

    #[inline(always)]
    pub fn neg(self) -> Fp {
        Fp::ZERO.sub(self)
    }

    #[inline(always)]
    pub fn mul(self, o: Fp) -> Fp {
        let prod = (self.0 as u128) * (o.0 as u128); // < 2^122
        let lo = (prod as u64) & P;
        let hi = (prod >> 61) as u64; // < 2^61
        let r = ct_sub_p(lo + hi); // lo + hi <= 2P
        Fp(ct_sub_p(r))
    }

    /// Constant-time select: returns `a` if choice == 1, `b` if choice == 0.
    #[inline(always)]
    pub fn select(choice: u64, a: Fp, b: Fp) -> Fp {
        let m = 0u64.wrapping_sub(choice & 1);
        Fp((a.0 & m) | (b.0 & !m))
    }

    /// Fixed-schedule exponentiation (61 squarings + 61 selects; no
    /// exponent- or base-dependent branches).
    pub fn pow(self, e: u64) -> Fp {
        let mut acc = Fp::ONE;
        let mut i = 64;
        while i > 0 {
            i -= 1;
            acc = acc.mul(acc);
            let m = acc.mul(self);
            acc = Fp::select((e >> i) & 1, m, acc);
        }
        acc
    }

    /// Multiplicative inverse via Fermat: a^(p-2). inv(0) = 0.
    pub fn inv(self) -> Fp {
        self.pow(P - 2)
    }

    /// Constant-time equality, returns 1 if equal else 0.
    #[inline(always)]
    pub fn ct_eq(self, o: Fp) -> u64 {
        let d = self.0 ^ o.0;
        // d == 0  ->  (d | -d) has top bit 0
        1 ^ ((d | d.wrapping_neg()) >> 63)
    }

    /// Variable-time zero test; use only on PUBLIC values.
    pub fn is_zero_public(self) -> bool {
        self.0 == 0
    }
}

impl PartialEq for Fp {
    fn eq(&self, o: &Fp) -> bool {
        self.ct_eq(*o) == 1
    }
}
impl Eq for Fp {}

/// Overwrite a slice of field elements with zero in a way the optimizer
/// may not elide.
pub fn zeroize(v: &mut [Fp]) {
    for x in v.iter_mut() {
        unsafe { core::ptr::write_volatile(&mut x.0, 0) };
    }
    compiler_fence(Ordering::SeqCst);
}

pub fn zeroize_bytes(v: &mut [u8]) {
    for x in v.iter_mut() {
        unsafe { core::ptr::write_volatile(x, 0) };
    }
    compiler_fence(Ordering::SeqCst);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_mul(a: u64, b: u64) -> u64 {
        ((a as u128 * b as u128) % (P as u128)) as u64
    }

    #[test]
    fn arithmetic_matches_reference() {
        let mut s: u64 = 0x9E3779B97F4A7C15;
        let edge = [0, 1, 2, P - 1, P - 2, 1 << 60, (1 << 60) + 1];
        let mut vals: Vec<u64> = edge.to_vec();
        for _ in 0..2000 {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            vals.push(s % P);
        }
        for i in 0..vals.len() {
            let a = vals[i];
            let b = vals[(i * 7 + 3) % vals.len()];
            let (fa, fb) = (Fp::new(a), Fp::new(b));
            assert_eq!(fa.mul(fb).0, naive_mul(a, b));
            assert_eq!(fa.add(fb).0, ((a as u128 + b as u128) % P as u128) as u64);
            assert_eq!(
                fa.sub(fb).0,
                ((a as u128 + P as u128 - b as u128) % P as u128) as u64
            );
            if a != 0 {
                assert_eq!(fa.mul(fa.inv()).0, 1);
            }
        }
        assert_eq!(Fp::new(u64::MAX).0, (u64::MAX % P));
        assert_eq!(Fp::new(P).0, 0);
    }

    #[test]
    fn ct_eq_works() {
        assert_eq!(Fp(5).ct_eq(Fp(5)), 1);
        assert_eq!(Fp(5).ct_eq(Fp(6)), 0);
        assert_eq!(Fp(0).ct_eq(Fp(P - 1)), 0);
    }
}
