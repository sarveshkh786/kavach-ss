//! Public Cauchy-matrix linear algebra.
//!
//! The extended generator is the (n+1) x t Cauchy matrix
//!     C[i][j] = 1 / (x_i - y_j),   x_i = i (i = 0..=n),   y_j = p - 1 - j (j = 0..t-1).
//! Every square submatrix of a Cauchy matrix is non-singular, which gives
//! both reconstruction (any t share-rows) and perfect privacy (row 0 plus
//! any t-1 share-rows). All data here is PUBLIC (indices only), so
//! variable-time Gaussian elimination is safe.

use crate::field::{Fp, P};

#[inline]
pub fn x(i: usize) -> Fp {
    Fp::new(i as u64)
}
#[inline]
pub fn y(j: usize) -> Fp {
    Fp::new(P - 1 - j as u64)
}

/// Generator entry C[i][j].
pub fn entry(i: usize, j: usize) -> Fp {
    x(i).sub(y(j)).inv()
}

/// Row i of the generator, length t.
pub fn row(i: usize, t: usize) -> Vec<Fp> {
    (0..t).map(|j| entry(i, j)).collect()
}

/// Solve A^T z = b for z where A is k x k (rows given). Returns None if singular.
/// Public data only.
fn solve_transposed(rows: &[Vec<Fp>], b: &[Fp]) -> Option<Vec<Fp>> {
    let k = rows.len();
    // Build augmented matrix M = [A^T | b]
    let mut m: Vec<Vec<Fp>> = (0..k)
        .map(|r| {
            let mut v: Vec<Fp> = (0..k).map(|c| rows[c][r]).collect();
            v.push(b[r]);
            v
        })
        .collect();
    for col in 0..k {
        let piv = (col..k).find(|&r| !m[r][col].is_zero_public())?;
        m.swap(col, piv);
        let inv = m[col][col].inv();
        for c in col..=k {
            m[col][c] = m[col][c].mul(inv);
        }
        for r in 0..k {
            if r != col && !m[r][col].is_zero_public() {
                let f = m[r][col];
                for c in col..=k {
                    let d = f.mul(m[col][c]);
                    m[r][c] = m[r][c].sub(d);
                }
            }
        }
    }
    Some(m.into_iter().map(|r| r[k]).collect())
}

/// Recombination coefficients: given t distinct share indices I (each in 1..=n),
/// returns mu such that  C_target . v  =  sum_k mu_k * share_{I_k}
/// for every v. (mu = C_target . C_I^{-1}).
pub fn recombination(indices: &[usize], target: usize, t: usize) -> Option<Vec<Fp>> {
    assert_eq!(indices.len(), t);
    let rows: Vec<Vec<Fp>> = indices.iter().map(|&i| row(i, t)).collect();
    let tgt = row(target, t);
    // mu . C_I = tgt  <=>  C_I^T mu^T = tgt^T
    solve_transposed(&rows, &tgt)
}

/// Rank test over F_p (public data). Used by tests and the `audit` command.
pub fn is_nonsingular(rows: &[Vec<Fp>]) -> bool {
    let k = rows.len();
    let mut m: Vec<Vec<Fp>> = rows.to_vec();
    for col in 0..k {
        let piv = match (col..k).find(|&r| !m[r][col].is_zero_public()) {
            Some(p) => p,
            None => return false,
        };
        m.swap(col, piv);
        let inv = m[col][col].inv();
        for r in (col + 1)..k {
            let f = m[r][col].mul(inv);
            for c in col..k {
                let d = f.mul(m[col][c]);
                m[r][c] = m[r][c].sub(d);
            }
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn subsets(n: usize, k: usize, start: usize, cur: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if cur.len() == k {
            out.push(cur.clone());
            return;
        }
        for i in start..=n {
            cur.push(i);
            subsets(n, k, i + 1, cur, out);
            cur.pop();
        }
    }

    /// The two linear-algebra facts the security proof rests on, checked
    /// exhaustively for all (t, n) with n <= 8.
    #[test]
    fn mds_and_privacy_conditions_hold() {
        for n in 1..=8usize {
            for t in 1..=n {
                // Reconstruction: every t-subset of share rows {1..n} is invertible.
                let mut subs = vec![];
                subsets(n, t, 1, &mut vec![], &mut subs);
                for s in &subs {
                    let rows: Vec<_> = s.iter().map(|&i| row(i, t)).collect();
                    assert!(is_nonsingular(&rows), "recon fails t={t} n={n} {s:?}");
                }
                // Privacy: row 0 plus any (t-1)-subset of share rows is invertible.
                let mut subs = vec![];
                subsets(n, t - 1, 1, &mut vec![], &mut subs);
                for s in &subs {
                    let mut rows = vec![row(0, t)];
                    rows.extend(s.iter().map(|&i| row(i, t)));
                    assert!(is_nonsingular(&rows), "privacy fails t={t} n={n} {s:?}");
                }
            }
        }
    }
}
