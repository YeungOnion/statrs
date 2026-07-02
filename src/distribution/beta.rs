use crate::function::{beta, gamma};
use crate::prec;

/// Implements the [Beta](https://en.wikipedia.org/wiki/Beta_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::{Beta, Continuous};
/// use statrs::statistics::*;
/// use approx::assert_abs_diff_eq;
///
/// let n = Beta::new(2.0, 2.0).unwrap();
/// assert_eq!(n.mean().unwrap(), 0.5);
/// assert_abs_diff_eq!(n.pdf(0.5), 1.5, epsilon = 1e-14);
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Beta {
    shape_a: f64,
    shape_b: f64,
}

/// Represents the errors that can occur when creating a [`Beta`].
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[non_exhaustive]
pub enum BetaError {
    /// Shape A is NaN, infinite, zero or negative.
    ShapeAInvalid,

    /// Shape B is NaN, infinite, zero or negative.
    ShapeBInvalid,
}

impl core::fmt::Display for BetaError {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            BetaError::ShapeAInvalid => write!(f, "Shape A is NaN, infinite, zero or negative"),
            BetaError::ShapeBInvalid => write!(f, "Shape B is NaN, infinite, zero or negative"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for BetaError {}

impl Beta {
    /// Constructs a new beta distribution with shapeA (α) of `shape_a`
    /// and shapeB (β) of `shape_b`
    ///
    /// # Errors
    ///
    /// Returns an error if `shape_a` or `shape_b` are `NaN` or infinite.
    /// Also returns an error if `shape_a <= 0.0` or `shape_b <= 0.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Beta;
    ///
    /// let mut result = Beta::new(2.0, 2.0);
    /// assert!(result.is_ok());
    ///
    /// result = Beta::new(0.0, 0.0);
    /// assert!(result.is_err());
    /// ```
    pub fn new(shape_a: f64, shape_b: f64) -> Result<Beta, BetaError> {
        if shape_a.is_nan() || shape_a.is_infinite() || shape_a <= 0.0 {
            return Err(BetaError::ShapeAInvalid);
        }

        if shape_b.is_nan() || shape_b.is_infinite() || shape_b <= 0.0 {
            return Err(BetaError::ShapeBInvalid);
        }

        Ok(Beta { shape_a, shape_b })
    }

    /// Returns the shapeA (α) of the beta distribution
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Beta;
    ///
    /// let n = Beta::new(1.0, 2.0).unwrap();
    /// assert_eq!(n.shape_a(), 1.0);
    /// ```
    pub fn shape_a(&self) -> f64 {
        self.shape_a
    }

    /// Returns the shapeB (β) of the beta distributionβ
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Beta;
    ///
    /// let n = Beta::new(1.0, 2.0).unwrap();
    /// assert_eq!(n.shape_b(), 2.0);
    /// ```
    pub fn shape_b(&self) -> f64 {
        self.shape_b
    }
}

mod legacy;

pub(crate) fn pdf(shape_a: f64, shape_b: f64, x: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        0.0
    } else if prec::ulps_eq!(shape_a, 1.0) && prec::ulps_eq!(shape_b, 1.0) {
        1.0
    } else if shape_a > 80.0 || shape_b > 80.0 {
        ln_pdf(shape_a, shape_b, x).exp()
    } else {
        let bb = gamma::gamma(shape_a + shape_b) / (gamma::gamma(shape_a) * gamma::gamma(shape_b));
        bb * x.powf(shape_a - 1.0) * (1.0 - x).powf(shape_b - 1.0)
    }
}

pub(crate) fn ln_pdf(shape_a: f64, shape_b: f64, x: f64) -> f64 {
    if !(0.0..=1.0).contains(&x) {
        f64::NEG_INFINITY
    } else if prec::ulps_eq!(shape_a, 1.0) && prec::ulps_eq!(shape_b, 1.0) {
        0.0
    } else {
        let aa = gamma::ln_gamma(shape_a + shape_b) - gamma::ln_gamma(shape_a) - gamma::ln_gamma(shape_b);
        let bb = if prec::ulps_eq!(shape_a, 1.0) && x == 0.0 {
            0.0
        } else if x == 0.0 {
            f64::NEG_INFINITY
        } else {
            (shape_a - 1.0) * x.ln()
        };
        let cc = if prec::ulps_eq!(shape_b, 1.0) && prec::ulps_eq!(x, 1.0) {
            0.0
        } else if prec::ulps_eq!(x, 1.0) {
            f64::NEG_INFINITY
        } else {
            (shape_b - 1.0) * (1.0 - x).ln()
        };
        aa + bb + cc
    }
}

pub(crate) fn cdf(shape_a: f64, shape_b: f64, x: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else if x >= 1.0 {
        1.0
    } else if prec::ulps_eq!(shape_a, 1.0) && prec::ulps_eq!(shape_b, 1.0) {
        x
    } else {
        beta::beta_reg(shape_a, shape_b, x)
    }
}

pub(crate) fn sf(shape_a: f64, shape_b: f64, x: f64) -> f64 {
    if x < 0.0 {
        1.0
    } else if x >= 1.0 {
        0.0
    } else if prec::ulps_eq!(shape_a, 1.0) && prec::ulps_eq!(shape_b, 1.0) {
        1. - x
    } else {
        beta::beta_reg(shape_b, shape_a, 1.0 - x)
    }
}

pub(crate) fn inverse_cdf_unchecked(shape_a: f64, shape_b: f64, p: f64) -> f64 {
    beta::inv_beta_reg(shape_a, shape_b, p)
}

pub(crate) fn mean(shape_a: f64, shape_b: f64) -> f64 {
    shape_a / (shape_a + shape_b)
}

pub(crate) fn variance(shape_a: f64, shape_b: f64) -> f64 {
    shape_a * shape_b / ((shape_a + shape_b) * (shape_a + shape_b) * (shape_a + shape_b + 1.0))
}

impl core::fmt::Display for Beta {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Beta(a={}, b={})", self.shape_a, self.shape_b)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distr::Distribution<f64> for Beta {
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        // Generated by sampling two gamma distributions and normalizing.
        let x = super::gamma::sample_unchecked(rng, self.shape_a, 1.0);
        let y = super::gamma::sample_unchecked(rng, self.shape_b, 1.0);
        x / (x + y)
    }
}

