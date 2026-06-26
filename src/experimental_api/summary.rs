use core::ops::ControlFlow;

use crate::experimental_api::distribution::Moments;
use crate::experimental_api::fold::FoldStat;
use crate::experimental_api::streaming::MeanAccum;

pub fn quadratic_mean(data: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for x in data {
        if x.is_nan() {
            return f64::NAN;
        }
        sum += x * x;
        count += 1;
    }
    if count == 0 {
        return f64::NAN;
    }
    (sum / count as f64).sqrt()
}

// --- geometric mean ---

pub fn geometric_step(acc: MeanAccum, x: f64) -> ControlFlow<Result<f64, f64>, MeanAccum> {
    if x > 0.0 {
        ControlFlow::Continue(acc.push(x.ln()))
    } else {
        ControlFlow::Break(Err(x))
    }
}

pub fn geometric_finalize(acc: MeanAccum) -> Option<f64> {
    acc.mean().map(f64::exp)
}

pub fn geometric_mean(iter: impl IntoIterator<Item = f64>) -> Result<Option<f64>, f64> {
    iter.into_iter()
        .fold_stat(MeanAccum::default(), geometric_step, geometric_finalize)
}

// --- harmonic mean ---

pub fn harmonic_step(acc: MeanAccum, x: f64) -> ControlFlow<Result<f64, f64>, MeanAccum> {
    if x > 0.0 {
        ControlFlow::Continue(acc.push(x.recip()))
    } else {
        ControlFlow::Break(Err(x))
    }
}

pub fn harmonic_finalize(acc: MeanAccum) -> Option<f64> {
    acc.mean().map(f64::recip)
}

pub fn harmonic_mean(iter: impl IntoIterator<Item = f64>) -> Result<Option<f64>, f64> {
    iter.into_iter()
        .fold_stat(MeanAccum::default(), harmonic_step, harmonic_finalize)
}

// --- abs_min / abs_max ---
//
// State is Option<f64>: None until the first element is seen.
// NaN is unordered and short-circuits as an error; all finite and
// infinite values are valid absolute-value comparisons.

pub fn abs_min_step(acc: Option<f64>, x: f64) -> ControlFlow<Result<f64, f64>, Option<f64>> {
    if x.is_nan() {
        ControlFlow::Break(Err(x))
    } else {
        let a = x.abs();
        ControlFlow::Continue(Some(acc.map_or(a, |m| m.min(a))))
    }
}

pub fn abs_min_finalize(acc: Option<f64>) -> Option<f64> {
    acc
}

pub fn abs_min(iter: impl IntoIterator<Item = f64>) -> Result<Option<f64>, f64> {
    iter.into_iter()
        .fold_stat(None, abs_min_step, abs_min_finalize)
}

pub fn abs_max_step(acc: Option<f64>, x: f64) -> ControlFlow<Result<f64, f64>, Option<f64>> {
    if x.is_nan() {
        ControlFlow::Break(Err(x))
    } else {
        let a = x.abs();
        ControlFlow::Continue(Some(acc.map_or(a, |m| m.max(a))))
    }
}

pub fn abs_max_finalize(acc: Option<f64>) -> Option<f64> {
    acc
}

pub fn abs_max(iter: impl IntoIterator<Item = f64>) -> Result<Option<f64>, f64> {
    iter.into_iter()
        .fold_stat(None, abs_max_step, abs_max_finalize)
}

pub mod steps {
    pub use super::{
        abs_max_finalize, abs_max_step, abs_min_finalize, abs_min_step, geometric_finalize,
        geometric_step, harmonic_finalize, harmonic_step,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    // quadratic_mean (unchanged behaviour)

    #[test]
    fn quadratic_mean_empty_is_nan() {
        assert!(quadratic_mean([] as [f64; 0]).is_nan());
    }

    #[test]
    fn quadratic_mean_nan_propagates() {
        assert!(quadratic_mean([1.0, f64::NAN, 2.0]).is_nan());
    }

    #[test]
    fn quadratic_mean_known_value() {
        let result = quadratic_mean([0.0, 3.0, -2.0]);
        assert!((result - 2.08167_f64).abs() < 1e-4);
    }

    // geometric_mean

    #[test]
    fn geometric_mean_empty() {
        assert_eq!(geometric_mean([]), Ok(None));
    }

    #[test]
    fn geometric_mean_known_value() {
        // [1, 2, 3]: geometric mean = (6)^(1/3) ≈ 1.81712
        let r = geometric_mean([1.0, 2.0, 3.0]).unwrap().unwrap();
        assert!((r - 1.81712_f64).abs() < 1e-4);
    }

    #[test]
    fn geometric_mean_non_positive_errors() {
        assert_eq!(geometric_mean([1.0, 0.0, 2.0]), Err(0.0));
        assert_eq!(geometric_mean([1.0, -1.0, 2.0]), Err(-1.0));
    }

    #[test]
    fn geometric_mean_nan_errors() {
        let r = geometric_mean([1.0, f64::NAN]);
        assert!(matches!(r, Err(x) if x.is_nan()));
    }

    // harmonic_mean

    #[test]
    fn harmonic_mean_empty() {
        assert_eq!(harmonic_mean([]), Ok(None));
    }

    #[test]
    fn harmonic_mean_known_value() {
        // [1, 2, 3]: harmonic mean = 3 / (1 + 1/2 + 1/3) = 18/11 ≈ 1.63636
        let r = harmonic_mean([1.0, 2.0, 3.0]).unwrap().unwrap();
        assert!((r - 18.0_f64 / 11.0).abs() < 1e-10);
    }

    #[test]
    fn harmonic_mean_non_positive_errors() {
        assert_eq!(harmonic_mean([1.0, 0.0, 2.0]), Err(0.0));
        assert_eq!(harmonic_mean([1.0, -1.0, 2.0]), Err(-1.0));
    }

    // abs_min

    #[test]
    fn abs_min_empty() {
        assert_eq!(abs_min([]), Ok(None));
    }

    #[test]
    fn abs_min_known_value() {
        assert_eq!(abs_min([3.0, -1.0, -4.0, 2.0]), Ok(Some(1.0)));
    }

    #[test]
    fn abs_min_nan_errors() {
        let r = abs_min([1.0, f64::NAN, 2.0]);
        assert!(matches!(r, Err(x) if x.is_nan()));
    }

    // abs_max

    #[test]
    fn abs_max_empty() {
        assert_eq!(abs_max([]), Ok(None));
    }

    #[test]
    fn abs_max_known_value() {
        assert_eq!(abs_max([3.0, -8.0, 2.0]), Ok(Some(8.0)));
    }

    #[test]
    fn abs_max_nan_errors() {
        let r = abs_max([1.0, f64::NAN, 2.0]);
        assert!(matches!(r, Err(x) if x.is_nan()));
    }
}
