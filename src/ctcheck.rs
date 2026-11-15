//! dudect-style leakage test (Reparaz, Balasch, Verbauwhede 2017):
//! time an operation on two input classes (fixed vs random secret) in
//! randomly interleaved order and apply Welch's t-test. |t| < 4.5 means no
//! timing difference was detected at that sample size.

use crate::field::Fp;
use crate::sharing::Rng;
use std::hint::black_box;
use std::time::Instant;

pub struct Report {
    pub name: &'static str,
    pub samples: usize,
    pub t: f64,
}

fn welch(a: &[f64], b: &[f64]) -> f64 {
    let m = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let var = |v: &[f64], mu: f64| {
        v.iter().map(|x| (x - mu).powi(2)).sum::<f64>() / (v.len() as f64 - 1.0)
    };
    let (ma, mb) = (m(a), m(b));
    let (va, vb) = (var(a, ma), var(b, mb));
    (ma - mb) / (va / a.len() as f64 + vb / b.len() as f64).sqrt()
}

/// Crop the top 10% of each class (OS-noise outliers) as dudect does.
fn crop(mut v: Vec<f64>) -> Vec<f64> {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let k = v.len() * 9 / 10;
    v.truncate(k);
    v
}

fn run<F: Fn(Fp, Fp) -> Fp>(name: &'static str, samples: usize, op: F, rng: &mut Rng) -> Report {
    const INNER: usize = 64;
    let mut a = Vec::with_capacity(samples);
    let mut b = Vec::with_capacity(samples);
    let fixed = Fp::new(0); // pathological class: zero secret
                            // Pre-generate every input so the measured region is identical for both
                            // classes (no RNG syscalls / cache pollution inside or right before it).
    let total = samples * 2;
    let mut cls = vec![0u8; total];
    rng.bytes(&mut cls).unwrap();
    let inputs: Vec<(Fp, Fp)> = (0..total)
        .map(|i| {
            let r = rng.fp().unwrap();
            let o = rng.fp().unwrap();
            (if cls[i] & 1 == 0 { fixed } else { r }, o)
        })
        .collect();
    // warm-up
    for &(s, o) in inputs.iter().take(1000) {
        black_box(op(black_box(s), black_box(o)));
    }
    for i in 0..total {
        let (secret, other) = inputs[i];
        let start = Instant::now();
        let mut acc = secret;
        for _ in 0..INNER {
            acc = black_box(op(black_box(acc), black_box(other)));
        }
        black_box(acc);
        let dt = start.elapsed().as_nanos() as f64;
        if cls[i] & 1 == 0 {
            a.push(dt)
        } else {
            b.push(dt)
        }
    }
    Report {
        name,
        samples: a.len() + b.len(),
        t: welch(&crop(a), &crop(b)),
    }
}

pub fn run_all(samples: usize) -> Vec<Report> {
    let mut rng = Rng::new().expect("rng");
    vec![
        run(
            "Fp::mul  (secret x public)",
            samples,
            |s, o| s.mul(o),
            &mut rng,
        ),
        run(
            "Fp::add  (secret + public)",
            samples,
            |s, o| s.add(o),
            &mut rng,
        ),
        run(
            "Fp::sub  (secret - public)",
            samples,
            |s, o| s.sub(o),
            &mut rng,
        ),
        run(
            "Fp::ct_eq (MAC tag compare)",
            samples,
            |s, o| Fp::new(s.ct_eq(o)).add(s),
            &mut rng,
        ),
        run(
            "Fp::pow  (fixed schedule)",
            samples / 8,
            |s, _o| s.pow(0x1F2E3D4C5B6A7988),
            &mut rng,
        ),
    ]
}
