use crate::consts;
use crate::function::erf;
use core::f64;

/// Implements the [Normal](https://en.wikipedia.org/wiki/Normal_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::Normal;
///
/// let n = Normal::new(0.0, 1.0).unwrap();
///
/// #[cfg(feature = "experimental_api")]
/// {
///     use statrs::experimental_api::{Mean, Pdf, TryVariate};
///
///     assert_eq!(n.mean(), 0.0);
///     let x = n.try_variate(1.0).unwrap();
///     assert_eq!(n.pdf(x).into_inner(), 0.2419707245191433497978);
/// }
///
/// #[cfg(not(feature = "experimental_api"))]
/// {
///     use statrs::distribution::Continuous;
///     use statrs::statistics::Distribution;
///
///     assert_eq!(n.mean().unwrap(), 0.0);
///     assert_eq!(n.pdf(1.0), 0.2419707245191433497978);
/// }
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Normal {
    mean: f64,
    std_dev: f64,
}

/// Represents the errors that can occur when creating a [`Normal`].
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[non_exhaustive]
pub enum NormalError {
    /// The mean is NaN.
    MeanInvalid,

    /// The standard deviation is NaN, zero or less than zero.
    StandardDeviationInvalid,
}

impl core::fmt::Display for NormalError {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            NormalError::MeanInvalid => write!(f, "Mean is NaN"),
            NormalError::StandardDeviationInvalid => {
                write!(f, "Standard deviation is NaN, zero or less than zero")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for NormalError {}

impl Normal {
    ///  Constructs a new normal distribution with a mean of `mean`
    /// and a standard deviation of `std_dev`
    ///
    /// # Errors
    ///
    /// Returns an error if `mean` or `std_dev` are `NaN` or if
    /// `std_dev <= 0.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Normal;
    ///
    /// let mut result = Normal::new(0.0, 1.0);
    /// assert!(result.is_ok());
    ///
    /// result = Normal::new(0.0, 0.0);
    /// assert!(result.is_err());
    /// ```
    pub fn new(mean: f64, std_dev: f64) -> Result<Normal, NormalError> {
        if mean.is_nan() {
            return Err(NormalError::MeanInvalid);
        }

        if std_dev.is_nan() || std_dev <= 0.0 {
            return Err(NormalError::StandardDeviationInvalid);
        }

        Ok(Normal { mean, std_dev })
    }

    /// Constructs a new standard normal distribution with a mean of 0
    /// and a standard deviation of 1.
    ///
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Normal;
    ///
    /// let mut result = Normal::standard();
    /// ```
    pub fn standard() -> Normal {
        Normal {
            mean: 0.0,
            std_dev: 1.0,
        }
    }
}

#[cfg(not(feature = "experimental_api"))]
mod legacy;

impl core::fmt::Display for Normal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "N({},{})", self.mean, self.std_dev)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distr::Distribution<f64> for Normal {
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        sample_unchecked(rng, self.mean, self.std_dev)
    }
}

/// performs an unchecked cdf calculation for a normal distribution
/// with the given mean and standard deviation at x
pub fn cdf_unchecked(x: f64, mean: f64, std_dev: f64) -> f64 {
    0.5 * erf::erfc((mean - x) / (std_dev * f64::consts::SQRT_2))
}

/// performs an unchecked sf calculation for a normal distribution
/// with the given mean and standard deviation at x
pub fn sf_unchecked(x: f64, mean: f64, std_dev: f64) -> f64 {
    0.5 * erf::erfc((x - mean) / (std_dev * f64::consts::SQRT_2))
}

/// performs an unchecked pdf calculation for a normal distribution
/// with the given mean and standard deviation at x
pub fn pdf_unchecked(x: f64, mean: f64, std_dev: f64) -> f64 {
    let d = (x - mean) / std_dev;
    (-0.5 * d * d).exp() / (consts::SQRT_2PI * std_dev)
}

/// performs an unchecked log(pdf) calculation for a normal distribution
/// with the given mean and standard deviation at x
pub fn ln_pdf_unchecked(x: f64, mean: f64, std_dev: f64) -> f64 {
    let d = (x - mean) / std_dev;
    (-0.5 * d * d) - consts::LN_SQRT_2PI - std_dev.ln()
}

/// performs an unchecked inverse cdf calculation for a normal distribution
/// with the given mean and standard deviation at p
pub fn inverse_cdf_unchecked(p: f64, mean: f64, std_dev: f64) -> f64 {
    mean - (std_dev * f64::consts::SQRT_2 * erf::erfc_inv(2.0 * p))
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
/// draws a sample from a normal distribution using the Box-Muller algorithm
pub fn sample_unchecked<R: ::rand::Rng + ?Sized>(rng: &mut R, mean: f64, std_dev: f64) -> f64 {
    use crate::distribution::ziggurat;

    mean + std_dev * ziggurat::sample_std_normal(rng)
}

impl core::default::Default for Normal {
    /// Returns the standard normal distribution with a mean of 0
    /// and a standard deviation of 1.
    fn default() -> Self {
        Self::standard()
    }
}

