//! Derived summary statistics for `f64` iterators.
//!
//! Each function makes a single pass. For multiple statistics in one pass,
//! compose [`Accumulate`] impls in a tuple and call [`Iterator::fold`] directly.
//!
//! ## NaN and missing data
//!
//! These functions let `NaN` short circuit to `Some(NaN)` returns as it is a poisoning operation.
//! The theme of this behavior is we won't knowingly _create_ a `NaN`.
//! So, missing data must be handled via another avenue, usually filtering for the Some<T: Float> gives the right semantic.
//! See function level docs for what `None` return means outside of empty data.

use core::f64;
use core::ops::ControlFlow;

use num_traits::Zero;

use crate::experimental_api::distribution::Mean;
use crate::experimental_api::streaming::{MeanAccum, RunningCov};

/// Convenience function for quadratic means.
///
/// If results are unexpected, read the NaN policy from this module before writing an issue.
pub fn quadratic_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    match iter
        .into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x.is_nan() {
                ControlFlow::Break(Some(f64::NAN))
            } else if x.is_infinite() {
                ControlFlow::Break(Some(f64::INFINITY))
            } else {
                ControlFlow::Continue(acc.push(x * x))
            }
        })
        .map_continue(|acc| acc.mean().map(f64::sqrt))
    {
        ControlFlow::Break(x) | ControlFlow::Continue(x) => x,
    }
}

/// Convenience function for geometric means.
///
/// If results are unexpected, read the NaN policy from this module before writing an issue.
pub fn geometric_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    match iter
        .into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x.is_nan() {
                ControlFlow::Break(Some(f64::NAN))
            } else if x == 0.0 || x.is_infinite() {
                ControlFlow::Break(Some(x))
            } else if x < 0.0 {
                ControlFlow::Break(None)
            } else {
                ControlFlow::Continue(acc.push(x.ln()))
            }
        })
        .map_continue(|acc| acc.mean().map(f64::exp))
    {
        ControlFlow::Break(x) | ControlFlow::Continue(x) => x,
    }
}

/// Convenience function for harmonic means.
///
/// If results are unexpected, read the NaN policy from this module before writing an issue.
pub fn harmonic_mean(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    match iter
        .into_iter()
        .try_fold(MeanAccum::default(), |acc, x| {
            if x.is_nan() {
                ControlFlow::Break(Some(f64::NAN))
            } else if x <= 0.0 || x.is_infinite() {
                ControlFlow::Break(None)
            } else {
                ControlFlow::Continue(acc.push(x.recip()))
            }
        })
        .map_continue(|acc| acc.mean().map(f64::recip))
    {
        ControlFlow::Break(x) | ControlFlow::Continue(x) => x,
    }
}

/// Convenience function for abs min.
///
/// If results are unexpected, read the NaN policy from this module before writing an issue.
pub fn abs_min(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    match iter.into_iter().try_fold(None::<f64>, |acc, x| {
        if x.is_nan() {
            ControlFlow::Break(Some(f64::NAN))
        } else if x.is_zero() {
            ControlFlow::Break(Some(0.0))
        } else {
            let a = x.abs();
            ControlFlow::Continue(Some(acc.map_or(a, |m| m.min(a))))
        }
    }) {
        ControlFlow::Break(x) | ControlFlow::Continue(x) => x,
    }
}

/// Convenience function for abs max.
///
/// If results are unexpected, read the NaN policy from this module before writing an issue.
pub fn abs_max(iter: impl IntoIterator<Item = f64>) -> Option<f64> {
    match iter.into_iter().try_fold(None::<f64>, |acc, x| {
        if x.is_nan() {
            ControlFlow::Break(Some(f64::NAN))
        } else if x.is_infinite() {
            ControlFlow::Break(Some(f64::INFINITY))
        } else {
            let a = x.abs();
            ControlFlow::Continue(Some(acc.map_or(a, |m| m.max(a))))
        }
    }) {
        ControlFlow::Break(x) | ControlFlow::Continue(x) => x,
    }
}

/// Sample covariance matrix over `N` variables from paired observations.
///
/// Returns `None` for fewer than two observations.
/// If results are unexpected, read the NaN policy from this module before writing an issue.
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
    fn quadratic_mean_nan_propagates() {
        assert!(quadratic_mean([1.0, f64::NAN, 2.0]).unwrap().is_nan());
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
    fn geometric_mean_zero_is_zero() {
        assert_eq!(geometric_mean([1.0, 0.0, 2.0]), Some(0.0));
    }

    #[test]
    fn geometric_mean_negative_is_none() {
        assert_eq!(geometric_mean([1.0, -1.0, 2.0]), None);
    }

    #[test]
    fn geometric_mean_nan_propagates() {
        assert!(geometric_mean([1.0, f64::NAN]).unwrap().is_nan());
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
    fn harmonic_mean_nan_propagates() {
        assert!(harmonic_mean([1.0, f64::NAN, 2.0]).unwrap().is_nan());
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
    fn abs_min_nan_propagates() {
        assert!(abs_min([1.0, f64::NAN, 2.0]).unwrap().is_nan());
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
    fn abs_max_nan_propagates() {
        assert!(abs_max([1.0, f64::NAN, 2.0]).unwrap().is_nan());
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
