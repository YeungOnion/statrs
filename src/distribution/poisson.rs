use crate::distribution::{Discrete, DiscreteCDF};
use crate::function::{factorial, gamma};
use crate::statistics::*;
use core::f64;

/// Implements the [Poisson](https://en.wikipedia.org/wiki/Poisson_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::{Poisson, Discrete};
/// use statrs::statistics::Distribution;
/// use approx::assert_abs_diff_eq;
///
/// let n = Poisson::new(1.0).unwrap();
/// assert_eq!(n.mean().unwrap(), 1.0);
/// assert_abs_diff_eq!(n.pmf(1), 0.367879441171442, epsilon = 1e-15);
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Poisson {
    lambda: f64,
}

/// Represents the errors that can occur when creating a [`Poisson`].
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[non_exhaustive]
pub enum PoissonError {
    /// The lambda is NaN, zero or less than zero.
    LambdaInvalid,
}

impl core::fmt::Display for PoissonError {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            PoissonError::LambdaInvalid => write!(f, "Lambda is NaN, zero or less than zero"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for PoissonError {}

impl Poisson {
    /// Constructs a new poisson distribution with a rate (λ)
    /// of `lambda`
    ///
    /// # Errors
    ///
    /// Returns an error if `lambda` is `NaN` or `lambda <= 0.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Poisson;
    ///
    /// let mut result = Poisson::new(1.0);
    /// assert!(result.is_ok());
    ///
    /// result = Poisson::new(0.0);
    /// assert!(result.is_err());
    /// ```
    pub fn new(lambda: f64) -> Result<Poisson, PoissonError> {
        if lambda.is_nan() || lambda <= 0.0 {
            Err(PoissonError::LambdaInvalid)
        } else {
            Ok(Poisson { lambda })
        }
    }

    /// Returns the rate (λ) of the poisson distribution
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Poisson;
    ///
    /// let n = Poisson::new(1.0).unwrap();
    /// assert_eq!(n.lambda(), 1.0);
    /// ```
    pub fn lambda(&self) -> f64 {
        self.lambda
    }
}

impl core::fmt::Display for Poisson {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Pois({})", self.lambda)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distributions::Distribution<u64> for Poisson {
    /// Generates one sample from the Poisson distribution either by
    /// Knuth's method if lambda < 30.0 or Rejection method PA by
    /// A. C. Atkinson from the Journal of the Royal Statistical Society
    /// Series C (Applied Statistics) Vol. 28 No. 1. (1979) pp. 29 - 35
    /// otherwise
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> u64 {
        sample_unchecked(rng, self.lambda) as u64
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distributions::Distribution<f64> for Poisson {
    /// Generates one sample from the Poisson distribution either by
    /// Knuth's method if lambda < 30.0 or Rejection method PA by
    /// A. C. Atkinson from the Journal of the Royal Statistical Society
    /// Series C (Applied Statistics) Vol. 28 No. 1. (1979) pp. 29 - 35
    /// otherwise
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        sample_unchecked(rng, self.lambda)
    }
}

impl DiscreteCDF<u64, f64> for Poisson {
    /// Calculates the cumulative distribution function for the poisson
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// P(x + 1, λ)
    /// ```
    ///
    /// where `λ` is the rate and `P` is the upper regularized gamma function
    fn cdf(&self, x: u64) -> f64 {
        gamma::gamma_ur(x as f64 + 1.0, self.lambda)
    }

    /// Calculates the survival function for the poisson
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// P(x + 1, λ)
    /// ```
    ///
    /// where `λ` is the rate and `P` is the lower regularized gamma function
    fn sf(&self, x: u64) -> f64 {
        gamma::gamma_lr(x as f64 + 1.0, self.lambda)
    }
}

impl Min<u64> for Poisson {
    /// Returns the minimum value in the domain of the poisson distribution
    /// representable by a 64-bit integer
    ///
    /// # Formula
    ///
    /// ```text
    /// 0
    /// ```
    fn min(&self) -> u64 {
        0
    }
}

impl Max<u64> for Poisson {
    /// Returns the maximum value in the domain of the poisson distribution
    /// representable by a 64-bit integer
    ///
    /// # Formula
    ///
    /// ```text
    /// 2^63 - 1
    /// ```
    fn max(&self) -> u64 {
        u64::MAX
    }
}

impl Distribution<f64> for Poisson {
    /// Returns the mean of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// λ
    /// ```
    ///
    /// where `λ` is the rate
    fn mean(&self) -> Option<f64> {
        Some(self.lambda)
    }

    /// Returns the variance of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// λ
    /// ```
    ///
    /// where `λ` is the rate
    fn variance(&self) -> Option<f64> {
        Some(self.lambda)
    }

    /// Returns the entropy of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// (1 / 2) * ln(2πeλ) - 1 / (12λ) - 1 / (24λ^2) - 19 / (360λ^3)
    /// ```
    ///
    /// where `λ` is the rate
    fn entropy(&self) -> Option<f64> {
        Some(
            0.5 * (2.0 * f64::consts::PI * f64::consts::E * self.lambda).ln()
                - 1.0 / (12.0 * self.lambda)
                - 1.0 / (24.0 * self.lambda * self.lambda)
                - 19.0 / (360.0 * self.lambda * self.lambda * self.lambda),
        )
    }

    /// Returns the skewness of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// λ^(-1/2)
    /// ```
    ///
    /// where `λ` is the rate
    fn skewness(&self) -> Option<f64> {
        Some(1.0 / self.lambda.sqrt())
    }
}

impl Median<f64> for Poisson {
    /// Returns the median of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// floor(λ + 1 / 3 - 0.02 / λ)
    /// ```
    ///
    /// where `λ` is the rate
    fn median(&self) -> f64 {
        (self.lambda + 1.0 / 3.0 - 0.02 / self.lambda).floor()
    }
}

impl Mode<Option<u64>> for Poisson {
    /// Returns the mode of the poisson distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// floor(λ)
    /// ```
    ///
    /// where `λ` is the rate
    fn mode(&self) -> Option<u64> {
        Some(self.lambda.floor() as u64)
    }
}

impl Discrete<u64, f64> for Poisson {
    /// Calculates the probability mass function for the poisson distribution at
    /// `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// (λ^x * e^(-λ)) / x!
    /// ```
    ///
    /// where `λ` is the rate
    fn pmf(&self, x: u64) -> f64 {
        (-self.lambda + x as f64 * self.lambda.ln() - factorial::ln_factorial(x)).exp()
    }

    /// Calculates the log probability mass function for the poisson
    /// distribution at
    /// `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// ln((λ^x * e^(-λ)) / x!)
    /// ```
    ///
    /// where `λ` is the rate
    fn ln_pmf(&self, x: u64) -> f64 {
        -self.lambda + x as f64 * self.lambda.ln() - factorial::ln_factorial(x)
    }
}

/// Generates one sample from the Poisson distribution either by
/// Knuth's method if lambda < 30.0 or Rejection method PA by
/// A. C. Atkinson from the Journal of the Royal Statistical Society
/// Series C (Applied Statistics) Vol. 28 No. 1. (1979) pp. 29 - 35
/// otherwise
#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
pub fn sample_unchecked<R: ::rand::Rng + ?Sized>(rng: &mut R, lambda: f64) -> f64 {
    if lambda < 30.0 {
        let limit = (-lambda).exp();
        let mut count = 0.0;
        let mut product: f64 = rng.gen();
        while product >= limit {
            count += 1.0;
            product *= rng.gen::<f64>();
        }
        count
    } else {
        let c = 0.767 - 3.36 / lambda;
        let beta = f64::consts::PI / (3.0 * lambda).sqrt();
        let alpha = beta * lambda;
        let k = c.ln() - lambda - beta.ln();

        loop {
            let u: f64 = rng.gen();
            let x = (alpha - ((1.0 - u) / u).ln()) / beta;
            let n = (x + 0.5).floor();
            if n < 0.0 {
                continue;
            }

            let v: f64 = rng.gen();
            let y = alpha - beta * x;
            let temp = 1.0 + y.exp();
            let lhs = y + (v / (temp * temp)).ln();
            let rhs = k + n * lambda.ln() - factorial::ln_factorial(n as u64);
            if lhs <= rhs {
                return n;
            }
        }
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::internal::density_util;
    use crate::distribution::internal::testing_boiler;
    testing_boiler!(lambda: f64; Poisson; PoissonError);

    #[test]
    fn test_create() {
        create_ok(1.5);
        create_ok(5.4);
        create_ok(10.8);
    }

    #[test]
    fn test_bad_create() {
        create_err(f64::NAN);
        create_err(-1.5);
        create_err(0.0);
    }

    #[test]
    fn test_mean() {
        let mean = |x: Poisson| x.mean().unwrap();
        test_exact(1.5, 1.5, mean);
        test_exact(5.4, 5.4, mean);
        test_exact(10.8, 10.8, mean);
    }

    #[test]
    fn test_variance() {
        let variance = |x: Poisson| x.variance().unwrap();
        test_exact(1.5, 1.5, variance);
        test_exact(5.4, 5.4, variance);
        test_exact(10.8, 10.8, variance);
    }

    #[test]
    fn test_entropy() {
        let entropy = |x: Poisson| x.entropy().unwrap();
        test_absolute(1.5, 1.531959153102376331946, 1e-15, entropy);
        test_absolute(5.4, 2.244941839577643504608, 1e-15, entropy);
        test_exact(10.8, 2.600596429676975222694, entropy);
    }

    #[test]
    fn test_skewness() {
        let skewness = |x: Poisson| x.skewness().unwrap();
        test_absolute(1.5, 0.8164965809277260327324, 1e-15, skewness);
        test_absolute(5.4, 0.4303314829119352094644, 1e-16, skewness);
        test_absolute(10.8, 0.3042903097250922852539, 1e-16, skewness);
    }

    #[test]
    fn test_median() {
        let median = |x: Poisson| x.median();
        test_exact(1.5, 1.0, median);
        test_exact(5.4, 5.0, median);
        test_exact(10.8, 11.0, median);
    }

    #[test]
    fn test_mode() {
        let mode = |x: Poisson| x.mode().unwrap();
        test_exact(1.5, 1, mode);
        test_exact(5.4, 5, mode);
        test_exact(10.8, 10, mode);
    }

    #[test]
    fn test_min_max() {
        let min = |x: Poisson| x.min();
        let max = |x: Poisson| x.max();
        test_exact(1.5, 0, min);
        test_exact(5.4, 0, min);
        test_exact(10.8, 0, min);
        test_exact(1.5, u64::MAX, max);
        test_exact(5.4, u64::MAX, max);
        test_exact(10.8, u64::MAX, max);
    }

    #[test]
    fn test_pmf() {
        let pmf = |arg: u64| move |x: Poisson| x.pmf(arg);
        test_absolute(1.5, 0.334695240222645000000000000000, 1e-15, pmf(1));
        test_absolute(1.5, 0.000003545747740570180000000000, 1e-20, pmf(10));
        test_absolute(1.5, 0.000000000000000304971208961018, 1e-30, pmf(20));
        test_absolute(5.4, 0.024389537090108400000000000000, 1e-17, pmf(1));
        test_absolute(5.4, 0.026241240591792300000000000000, 1e-16, pmf(10));
        test_absolute(5.4, 0.000000825202200316548000000000, 1e-20, pmf(20));
        test_absolute(10.8, 0.000220314636840657000000000000, 1e-18, pmf(1));
        test_absolute(10.8, 0.121365183659420000000000000000, 1e-15, pmf(10));
        test_absolute(10.8, 0.003908139778574110000000000000, 1e-16, pmf(20));
    }

    #[test]
    fn test_ln_pmf() {
        let ln_pmf = |arg: u64| move |x: Poisson| x.ln_pmf(arg);
        test_absolute(1.5, -1.09453489189183485135413967177, 1e-15, ln_pmf(1));
        test_absolute(1.5, -12.5497614919938728510400000000, 1e-14, ln_pmf(10));
        test_absolute(1.5, -35.7263142985901000000000000000, 1e-13, ln_pmf(20));
        test_exact(5.4, -3.71360104642977159156055355910, ln_pmf(1));
        test_absolute(5.4, -3.64042303737322774736223038530, 1e-15, ln_pmf(10));
        test_absolute(5.4, -14.0076373893489089949388000000, 1e-14, ln_pmf(20));
        test_absolute(10.8, -8.42045386586982559781714423000, 1e-14, ln_pmf(1));
        test_absolute(10.8, -2.10895123177378079525424989992, 1e-14, ln_pmf(10));
        test_absolute(10.8, -5.54469377815000936289610059500, 1e-14, ln_pmf(20));
    }

    #[test]
    fn test_cdf() {
        let cdf = |arg: u64| move |x: Poisson| x.cdf(arg);
        test_absolute(1.5, 0.5578254003710750000000, 1e-15, cdf(1));
        test_absolute(1.5, 0.9999994482467640000000, 1e-15, cdf(10));
        test_exact(1.5, 1.0, cdf(20));
        test_absolute(5.4, 0.0289061180327211000000, 1e-16, cdf(1));
        test_absolute(5.4, 0.9774863006897650000000, 1e-15, cdf(10));
        test_absolute(5.4, 0.9999997199928290000000, 1e-15, cdf(20));
        test_absolute(10.8, 0.0002407141402518290000, 1e-16, cdf(1));
        test_absolute(10.8, 0.4839692359955690000000, 1e-15, cdf(10));
        test_absolute(10.8, 0.9961800769608090000000, 1e-15, cdf(20));
    }

    #[test]
    fn test_sf() {
        let sf = |arg: u64| move |x: Poisson| x.sf(arg);
        test_absolute(1.5, 0.44217459962892536, 1e-15, sf(1));
        test_absolute(1.5, 0.0000005517532358246565, 1e-15, sf(10));
        test_absolute(1.5, 2.3372210700347092e-17, 1e-15, sf(20));
        test_absolute(5.4, 0.971093881967279, 1e-16, sf(1));
        test_absolute(5.4, 0.022513699310235582, 1e-15, sf(10));
        test_absolute(5.4, 0.0000002800071708975261, 1e-15, sf(20));
        test_absolute(10.8, 0.9997592858597482, 1e-16, sf(1));
        test_absolute(10.8, 0.5160307640044303, 1e-15, sf(10));
        test_absolute(10.8, 0.003819923039191422, 1e-15, sf(20));
    }

    #[test]
    fn test_discrete() {
        density_util::check_discrete_distribution(&create_ok(0.3), 10);
        density_util::check_discrete_distribution(&create_ok(4.5), 30);
    }

    #[test]
    fn test_inverse_cdf() {
        let invcdf = |arg: f64| move |x: Poisson| x.inverse_cdf(arg);
        test_exact(1.5, 0, invcdf(0.));
    }
}
