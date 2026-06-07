//! CNV caller — the prediction side of PyPGx's `Model[CNV]`
//! (`OneVsRestClassifier(SVC(kernel='rbf'))`).
//!
//! PyPGx ships these models as **pickled sklearn objects** that cannot load in
//! Rust. Path (A): extract each binary estimator's fitted params (support
//! vectors, dual coefficients, intercept, `gamma`) once in Python, then evaluate
//! the RBF decision function natively here — which reproduces sklearn's
//! `predict` **exactly** (verified in `tests/cnv.rs` against a reference model).
//! Training (`train_cnv_caller`) lives behind the `cnv` feature (needs sklears).

/// `scipy.ndimage.median_filter(c, size, mode='reflect')` for a 1-D signal.
/// Each output `i` is the rank-`size/2` element (lower-median for even `size`)
/// of the reflect-padded window `c[i - size/2 ..= i + ceil(size/2) - 1]`.
pub fn median_filter(c: &[f64], size: usize) -> Vec<f64> {
    let n = c.len();
    if n == 0 || size == 0 {
        return c.to_vec();
    }
    let half = size / 2;
    let reflect = |mut j: isize| -> usize {
        // scipy 'reflect' (half-sample symmetric): -1→0, n→n-1, …
        loop {
            if j < 0 {
                j = -j - 1;
            } else if j >= n as isize {
                j = 2 * n as isize - j - 1;
            } else {
                return j as usize;
            }
        }
    };
    // Reflect-pad once: `padded[k]` covers original index `k - half`, so output
    // `i` is the rank-`half` element of `padded[i .. i + size]`.
    let padded: Vec<f64> = (0..n + size - 1)
        .map(|k| c[reflect(k as isize - half as isize)])
        .collect();
    // Slide a kept-sorted window across `padded` (no per-window re-sort).
    let mut window: Vec<f64> = padded[..size].to_vec();
    window.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut out = Vec::with_capacity(n);
    out.push(window[half]);
    for i in 1..n {
        let rem = padded[i - 1];
        let pos = window
            .binary_search_by(|x| x.partial_cmp(&rem).unwrap())
            .expect("outgoing value present in window");
        window.remove(pos);
        let add = padded[i + size - 1];
        let ins = window
            .binary_search_by(|x| x.partial_cmp(&add).unwrap())
            .unwrap_or_else(|e| e);
        window.insert(ins, add);
        out.push(window[half]);
    }
    out
}

use serde::{Deserialize, Serialize};

/// One binary RBF-SVM estimator (`OneVsRest`: this class vs the rest).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CnvEstimator {
    pub label: i64,
    pub gamma: f64,
    pub intercept: f64,
    pub dual_coef: Vec<f64>,
    pub support_vectors: Vec<Vec<f64>>,
}

impl CnvEstimator {
    /// Binary SVC `decision_function(x)` = Σ αᵢ·K(svᵢ, x) + b, RBF kernel
    /// `K(u,v) = exp(-γ·‖u−v‖²)`.
    pub fn decision(&self, x: &[f64]) -> f64 {
        let mut acc = self.intercept;
        for (coef, sv) in self.dual_coef.iter().zip(&self.support_vectors) {
            let sq: f64 = sv.iter().zip(x).map(|(a, b)| (a - b) * (a - b)).sum();
            acc += coef * (-self.gamma * sq).exp();
        }
        acc
    }
}

/// A fitted `OneVsRestClassifier(SVC)` model. `estimators` are in `classes`
/// order; `predict` returns the class label whose estimator scores highest.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CnvModel {
    pub classes: Vec<i64>,
    pub estimators: Vec<CnvEstimator>,
}

impl CnvModel {
    /// `OneVsRestClassifier.predict(x)` — argmax of the per-class decision
    /// functions, returning the winning class label.
    pub fn predict(&self, x: &[f64]) -> i64 {
        let mut best = 0usize;
        let mut best_score = f64::NEG_INFINITY;
        for (k, est) in self.estimators.iter().enumerate() {
            let s = est.decision(x);
            if s > best_score {
                best_score = s;
                best = k;
            }
        }
        self.classes[best]
    }
}
