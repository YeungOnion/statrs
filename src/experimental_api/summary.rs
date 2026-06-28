//! Derived summary statistics for `f64` iterators.
//!
//! Each function makes a single pass. For multiple statistics in one pass,
//! compose [`Accumulate`] impls in a tuple and call [`Iterator::fold`] directly.
//!
//! [`Accumulate`]: crate::experimental_api::Accumulate

use crate::experimental_api::distribution::Moments;
use crate::experimental_api::streaming::{MeanAccum, RunningCov};

pub fn quadratic_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    iter.into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x.is_nan() { Err(()) } else { Ok(acc.push(x * x)) }
        })
        .ok()
        .and_then(|acc| acc.mean())
        .map(f64::sqrt)
}

pub fn geometric_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    iter.into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x > 0.0 { Ok(acc.push(x.ln())) } else { Err(()) }
        })
        .ok()
        .and_then(|acc| acc.mean())
        .map(f64::exp)
}

pub fn harmonic_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    iter.into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x > 0.0 { Ok(acc.push(x.recip())) } else { Err(()) }
        })
        .ok()
        .and_then(|acc| acc.mean())
        .map(f64::recip)
}

pub fn abs_min(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    iter.into_iter()
        .try_fold(None::<f64>, |acc, x| {
            if x.is_nan() {
                Err(())
            } else {
                let a = x.abs();
                Ok(Some(acc.map_or(a, |m| m.min(a))))
            }
        })
        .ok()
        .flatten()
}

pub fn abs_max(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    iter.into_iter()
        .try_fold(None::<f64>, |acc, x| {
            if x.is_nan() {
                Err(())
            } else {
                let a = x.abs();
                Ok(Some(acc.map_or(a, |m| m.max(a))))
            }
        })
        .ok()
        .flatten()
}

/// Sample covariance matrix over `N` variables from paired observations.
///
/// Returns `None` for fewer than two observations.
pub fn covariance<const N: usize>(
    data: impl IntoIterator<Item = [f64; N]>,
) -> Option<[[f64; N]; N]> {
    data.into_iter()
        .fold(RunningCov::<N>::default(), RunningCov::push)
        .finalize()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quadratic_mean_empty() {
        assert_eq!(quadratic_mean([] as [f64; 0]), None);
    }

    #[test]
    fn quadratic_mean_nan_is_none() {
        assert_eq!(quadratic_mean([1.0, f64::NAN, 2.0]), None);
    }

    #[test]
    fn quadratic_mean_known_value() {
        let r = quadratic_mean([0.0, 3.0, -2.0]).unwrap();
        assert!((r - 2.08167_f64).abs() < 1e-4);
    }

    #[test]
    fn geometric_mean_empty() {
        assert_eq!(geometric_mean([]), None);
    }

    #[test]
    fn geometric_mean_known_value() {
        let r = geometric_mean([1.0, 2.0, 3.0]).unwrap();
        assert!((r - 1.81712_f64).abs() < 1e-4);
    }

    #[test]
    fn geometric_mean_non_positive_is_none() {
        assert_eq!(geometric_mean([1.0, 0.0, 2.0]), None);
        assert_eq!(geometric_mean([1.0, -1.0, 2.0]), None);
    }

    #[test]
    fn geometric_mean_nan_is_none() {
        assert_eq!(geometric_mean([1.0, f64::NAN]), None);
    }

    #[test]
    fn harmonic_mean_empty() {
        assert_eq!(harmonic_mean([]), None);
    }

    #[test]
    fn harmonic_mean_known_value() {
        let r = harmonic_mean([1.0, 2.0, 3.0]).unwrap();
        assert!((r - 18.0_f64 / 11.0).abs() < 1e-10);
    }

    #[test]
    fn harmonic_mean_non_positive_is_none() {
        assert_eq!(harmonic_mean([1.0, 0.0, 2.0]), None);
        assert_eq!(harmonic_mean([1.0, -1.0, 2.0]), None);
    }

    #[test]
    fn abs_min_empty() {
        assert_eq!(abs_min([]), None);
    }

    #[test]
    fn abs_min_known_value() {
        assert_eq!(abs_min([3.0, -1.0, -4.0, 2.0]), Some(1.0));
    }

    #[test]
    fn abs_min_nan_is_none() {
        assert_eq!(abs_min([1.0, f64::NAN, 2.0]), None);
    }

    #[test]
    fn abs_max_empty() {
        assert_eq!(abs_max([]), None);
    }

    #[test]
    fn abs_max_known_value() {
        assert_eq!(abs_max([3.0, -8.0, 2.0]), Some(8.0));
    }

    #[test]
    fn abs_max_nan_is_none() {
        assert_eq!(abs_max([1.0, f64::NAN, 2.0]), None);
    }

    #[test]
    fn covariance_empty_is_none() {
        assert!(covariance::<2>([]).is_none());
    }

    #[test]
    fn covariance_n2_positive_correlation() {
        let data = [[1.0_f64, 2.0], [3.0, 4.0], [5.0, 6.0]];
        let mat = covariance(data).unwrap();
        assert!((mat[0][1] - 4.0).abs() < 1e-12);
    }
}
