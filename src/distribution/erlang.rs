use crate::distribution::{Continuous, ContinuousCDF, Gamma, GammaError};
use crate::statistics::*;

/// Implements the [Erlang](https://en.wikipedia.org/wiki/Erlang_distribution)
/// distribution
/// which is a special case of the
/// [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::{Erlang, Continuous};
/// use statrs::statistics::Distribution;
/// use approx::assert_abs_diff_eq;
///
/// let n = Erlang::new(3, 1.0).unwrap();
/// assert_eq!(n.mean().unwrap(), 3.0);
/// assert_abs_diff_eq!(n.pdf(2.0), 0.270670566473225383788, epsilon = 1e-15);
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Erlang {
    g: Gamma,
}

impl Erlang {
    /// Constructs a new erlang distribution with a shape (k)
    /// of `shape` and a rate (λ) of `rate`
    ///
    /// # Errors
    ///
    /// Returns an error if `shape` or `rate` are `NaN`.
    /// Also returns an error if `shape == 0` or `rate <= 0.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Erlang;
    ///
    /// let mut result = Erlang::new(3, 1.0);
    /// assert!(result.is_ok());
    ///
    /// result = Erlang::new(0, 0.0);
    /// assert!(result.is_err());
    /// ```
    pub fn new(shape: u64, rate: f64) -> Result<Erlang, GammaError> {
        Gamma::new(shape as f64, rate).map(|g| Erlang { g })
    }

    /// Returns the shape (k) of the erlang distribution
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Erlang;
    ///
    /// let n = Erlang::new(3, 1.0).unwrap();
    /// assert_eq!(n.shape(), 3);
    /// ```
    pub fn shape(&self) -> u64 {
        self.g.shape() as u64
    }

    /// Returns the rate (λ) of the erlang distribution
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Erlang;
    ///
    /// let n = Erlang::new(3, 1.0).unwrap();
    /// assert_eq!(n.rate(), 1.0);
    /// ```
    pub fn rate(&self) -> f64 {
        self.g.rate()
    }
}

impl core::fmt::Display for Erlang {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "E({}, {})", self.rate(), self.shape())
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distributions::Distribution<f64> for Erlang {
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        ::rand::distributions::Distribution::sample(&self.g, rng)
    }
}

impl ContinuousCDF<f64, f64> for Erlang {
    /// Calculates the cumulative distribution function for the erlang
    /// distribution
    /// at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// γ(k, λx)  (k - 1)!
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, and `γ` is the lower
    /// incomplete gamma function
    fn cdf(&self, x: f64) -> f64 {
        self.g.cdf(x)
    }

    /// Calculates the cumulative distribution function for the erlang
    /// distribution
    /// at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// γ(k, λx)  (k - 1)!
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, and `γ` is the upper
    /// incomplete gamma function
    fn sf(&self, x: f64) -> f64 {
        self.g.sf(x)
    }

    /// Calculates the inverse cumulative distribution function for the erlang
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// γ^{-1}(k, (k - 1)! x) / λ
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, and `γ` is the upper
    /// incomplete gamma function
    fn inverse_cdf(&self, p: f64) -> f64 {
        self.g.inverse_cdf(p)
    }
}

impl Min<f64> for Erlang {
    /// Returns the minimum value in the domain of the
    /// erlang distribution representable by a double precision
    /// float
    ///
    /// # Formula
    ///
    /// ```text
    /// 0
    /// ```
    fn min(&self) -> f64 {
        self.g.min()
    }
}

impl Max<f64> for Erlang {
    /// Returns the maximum value in the domain of the
    /// erlang distribution representable by a double precision
    /// float
    ///
    /// # Formula
    ///
    /// ```text
    /// f64::INFINITY
    /// ```
    fn max(&self) -> f64 {
        self.g.max()
    }
}

impl Distribution<f64> for Erlang {
    /// Returns the mean of the erlang distribution
    ///
    /// # Remarks
    ///
    /// Returns `shape` if `rate == f64::INFINITY`. This behavior
    /// is borrowed from the Math.NET implementation
    ///
    /// # Formula
    ///
    /// ```text
    /// k / λ
    /// ```
    ///
    /// where `k` is the shape and `λ` is the rate
    fn mean(&self) -> Option<f64> {
        self.g.mean()
    }

    /// Returns the variance of the erlang distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// k / λ^2
    /// ```
    ///
    /// where `α` is the shape and `λ` is the rate
    fn variance(&self) -> Option<f64> {
        self.g.variance()
    }

    /// Returns the entropy of the erlang distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// k - ln(λ) + ln(Γ(k)) + (1 - k) * ψ(k)
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, `Γ` is the gamma function,
    /// and `ψ` is the digamma function
    fn entropy(&self) -> Option<f64> {
        self.g.entropy()
    }

    /// Returns the skewness of the erlang distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// 2 / sqrt(k)
    /// ```
    ///
    /// where `k` is the shape
    fn skewness(&self) -> Option<f64> {
        self.g.skewness()
    }
}

impl Mode<Option<f64>> for Erlang {
    /// Returns the mode for the erlang distribution
    ///
    /// # Remarks
    ///
    /// Returns `shape` if `rate ==f64::INFINITY`. This behavior
    /// is borrowed from the Math.NET implementation
    ///
    /// # Formula
    ///
    /// ```text
    /// (k - 1) / λ
    /// ```
    ///
    /// where `k` is the shape and `λ` is the rate
    fn mode(&self) -> Option<f64> {
        self.g.mode()
    }
}

impl Continuous<f64, f64> for Erlang {
    /// Calculates the probability density function for the erlang distribution
    /// at `x`
    ///
    /// # Remarks
    ///
    /// Returns `NAN` if any of `shape` or `rate` are `INF`
    /// or if `x` is `INF`
    ///
    /// # Formula
    ///
    /// ```text
    /// (λ^k / Γ(k)) * x^(k - 1) * e^(-λ * x)
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, and `Γ` is the gamma function
    fn pdf(&self, x: f64) -> f64 {
        self.g.pdf(x)
    }

    /// Calculates the log probability density function for the erlang
    /// distribution
    /// at `x`
    ///
    /// # Remarks
    ///
    /// Returns `NAN` if any of `shape` or `rate` are `INF`
    /// or if `x` is `INF`
    ///
    /// # Formula
    ///
    /// ```text
    /// ln((λ^k / Γ(k)) * x^(k - 1) * e ^(-λ * x))
    /// ```
    ///
    /// where `k` is the shape, `λ` is the rate, and `Γ` is the gamma function
    fn ln_pdf(&self, x: f64) -> f64 {
        self.g.ln_pdf(x)
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::internal::density_util;
    use crate::distribution::internal::testing_boiler;
    testing_boiler!(shape: u64, rate: f64; Erlang; GammaError);

    #[test]
    fn test_create() {
        create_ok(1, 0.1);
        create_ok(1, 1.0);
        create_ok(10, 10.0);
        create_ok(10, 1.0);
        create_ok(10, f64::INFINITY);
    }

    #[test]
    fn test_bad_create() {
        let invalid = [
            (0, 1.0, GammaError::ShapeInvalid),
            (1, 0.0, GammaError::RateInvalid),
            (1, f64::NAN, GammaError::RateInvalid),
            (1, -1.0, GammaError::RateInvalid),
        ];

        for (s, r, err) in invalid {
            test_create_err(s, r, err);
        }
    }

    #[test]
    fn test_continuous() {
        density_util::check_continuous_distribution(&create_ok(1, 2.5), 0.0, 20.0);
        density_util::check_continuous_distribution(&create_ok(2, 1.5), 0.0, 20.0);
        density_util::check_continuous_distribution(&create_ok(3, 0.5), 0.0, 20.0);
    }
}
