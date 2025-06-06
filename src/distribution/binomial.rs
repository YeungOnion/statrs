use crate::distribution::{Discrete, DiscreteCDF};
use crate::function::{beta, factorial};
use crate::prec;
use crate::statistics::*;
use core::f64;

/// Implements the
/// [Binomial](https://en.wikipedia.org/wiki/Binomial_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::{Binomial, Discrete};
/// use statrs::statistics::Distribution;
///
/// let n = Binomial::new(0.5, 5).unwrap();
/// assert_eq!(n.mean().unwrap(), 2.5);
/// assert_eq!(n.pmf(0), 0.03125);
/// assert_eq!(n.pmf(3), 0.3125);
/// ```
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct Binomial {
    p: f64,
    n: u64,
}

/// Represents the errors that can occur when creating a [`Binomial`].
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[non_exhaustive]
pub enum BinomialError {
    /// The probability is NaN or not in `[0, 1]`.
    ProbabilityInvalid,
}

impl core::fmt::Display for BinomialError {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            BinomialError::ProbabilityInvalid => write!(f, "Probability is NaN or not in [0, 1]"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for BinomialError {}

impl Binomial {
    /// Constructs a new binomial distribution
    /// with a given `p` probability of success of `n`
    /// trials.
    ///
    /// # Errors
    ///
    /// Returns an error if `p` is `NaN`, less than `0.0`,
    /// greater than `1.0`, or if `n` is less than `0`
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Binomial;
    ///
    /// let mut result = Binomial::new(0.5, 5);
    /// assert!(result.is_ok());
    ///
    /// result = Binomial::new(-0.5, 5);
    /// assert!(result.is_err());
    /// ```
    pub fn new(p: f64, n: u64) -> Result<Binomial, BinomialError> {
        if p.is_nan() || !(0.0..=1.0).contains(&p) {
            Err(BinomialError::ProbabilityInvalid)
        } else {
            Ok(Binomial { p, n })
        }
    }

    /// Returns the probability of success `p` of
    /// the binomial distribution.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Binomial;
    ///
    /// let n = Binomial::new(0.5, 5).unwrap();
    /// assert_eq!(n.p(), 0.5);
    /// ```
    pub fn p(&self) -> f64 {
        self.p
    }

    /// Returns the number of trials `n` of the
    /// binomial distribution.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::distribution::Binomial;
    ///
    /// let n = Binomial::new(0.5, 5).unwrap();
    /// assert_eq!(n.n(), 5);
    /// ```
    pub fn n(&self) -> u64 {
        self.n
    }
}

impl core::fmt::Display for Binomial {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Bin({},{})", self.p, self.n)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distributions::Distribution<u64> for Binomial {
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> u64 {
        (0..self.n).fold(0, |acc, _| {
            let n: f64 = rng.gen();
            if n < self.p {
                acc + 1
            } else {
                acc
            }
        })
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distributions::Distribution<f64> for Binomial {
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        rng.sample::<u64, _>(self) as f64
    }
}

impl DiscreteCDF<u64, f64> for Binomial {
    /// Calculates the cumulative distribution function for the
    /// binomial distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// I_(1 - p)(n - x, 1 + x)
    /// ```
    ///
    /// where `I_(x)(a, b)` is the regularized incomplete beta function
    fn cdf(&self, x: u64) -> f64 {
        if x >= self.n {
            1.0
        } else {
            let k = x;
            beta::beta_reg((self.n - k) as f64, k as f64 + 1.0, 1.0 - self.p)
        }
    }

    /// Calculates the survival function for the
    /// binomial distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// I_(p)(x + 1, n - x)
    /// ```
    ///
    /// where `I_(x)(a, b)` is the regularized incomplete beta function
    fn sf(&self, x: u64) -> f64 {
        if x >= self.n {
            0.0
        } else {
            let k = x;
            beta::beta_reg(k as f64 + 1.0, (self.n - k) as f64, self.p)
        }
    }
}

impl Min<u64> for Binomial {
    /// Returns the minimum value in the domain of the
    /// binomial distribution representable by a 64-bit
    /// integer
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

impl Max<u64> for Binomial {
    /// Returns the maximum value in the domain of the
    /// binomial distribution representable by a 64-bit
    /// integer
    ///
    /// # Formula
    ///
    /// ```text
    /// n
    /// ```
    fn max(&self) -> u64 {
        self.n
    }
}

impl Distribution<f64> for Binomial {
    /// Returns the mean of the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// p * n
    /// ```
    fn mean(&self) -> Option<f64> {
        Some(self.p * self.n as f64)
    }

    /// Returns the variance of the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// n * p * (1 - p)
    /// ```
    fn variance(&self) -> Option<f64> {
        Some(self.p * (1.0 - self.p) * self.n as f64)
    }

    /// Returns the entropy of the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// (1 / 2) * ln (2 * π * e * n * p * (1 - p))
    /// ```
    fn entropy(&self) -> Option<f64> {
        let entr = if self.p == 0.0 || prec::ulps_eq!(self.p, 1.0) {
            0.0
        } else {
            (0..self.n + 1).fold(0.0, |acc, x| {
                let p = self.pmf(x);
                acc - p * p.ln()
            })
        };
        Some(entr)
    }

    /// Returns the skewness of the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// (1 - 2p) / sqrt(n * p * (1 - p)))
    /// ```
    fn skewness(&self) -> Option<f64> {
        Some((1.0 - 2.0 * self.p) / (self.n as f64 * self.p * (1.0 - self.p)).sqrt())
    }
}

impl Median<f64> for Binomial {
    /// Returns the median of the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// floor(n * p)
    /// ```
    fn median(&self) -> f64 {
        (self.p * self.n as f64).floor()
    }
}

impl Mode<Option<u64>> for Binomial {
    /// Returns the mode for the binomial distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// floor((n + 1) * p)
    /// ```
    fn mode(&self) -> Option<u64> {
        let mode = if self.p == 0.0 {
            0
        } else if prec::ulps_eq!(self.p, 1.0) {
            self.n
        } else {
            ((self.n as f64 + 1.0) * self.p).floor() as u64
        };
        Some(mode)
    }
}

impl Discrete<u64, f64> for Binomial {
    /// Calculates the probability mass function for the binomial
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// (n choose k) * p^k * (1 - p)^(n - k)
    /// ```
    fn pmf(&self, x: u64) -> f64 {
        if x > self.n {
            0.0
        } else if self.p == 0.0 {
            if x == 0 {
                1.0
            } else {
                0.0
            }
        } else if prec::ulps_eq!(self.p, 1.0) {
            if x == self.n {
                1.0
            } else {
                0.0
            }
        } else {
            (factorial::ln_binomial(self.n, x)
                + x as f64 * self.p.ln()
                + (self.n - x) as f64 * (1.0 - self.p).ln())
            .exp()
        }
    }

    /// Calculates the log probability mass function for the binomial
    /// distribution at `x`
    ///
    /// # Formula
    ///
    /// ```text
    /// ln((n choose k) * p^k * (1 - p)^(n - k))
    /// ```
    fn ln_pmf(&self, x: u64) -> f64 {
        if x > self.n {
            f64::NEG_INFINITY
        } else if self.p == 0.0 {
            if x == 0 {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        } else if prec::ulps_eq!(self.p, 1.0) {
            if x == self.n {
                0.0
            } else {
                f64::NEG_INFINITY
            }
        } else {
            factorial::ln_binomial(self.n, x)
                + x as f64 * self.p.ln()
                + (self.n - x) as f64 * (1.0 - self.p).ln()
        }
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::internal::density_util;
    crate::distribution::internal::testing_boiler!(p: f64, n: u64; Binomial; BinomialError);

    #[test]
    fn test_create() {
        create_ok(0.0, 4);
        create_ok(0.3, 3);
        create_ok(1.0, 2);
    }

    #[test]
    fn test_bad_create() {
        create_err(f64::NAN, 1);
        create_err(-1.0, 1);
        create_err(2.0, 1);
    }

    #[test]
    fn test_mean() {
        let mean = |x: Binomial| x.mean().unwrap();
        test_exact(0.0, 4, 0.0, mean);
        test_absolute(0.3, 3, 0.9, 1e-15, mean);
        test_exact(1.0, 2, 2.0, mean);
    }

    #[test]
    fn test_variance() {
        let variance = |x: Binomial| x.variance().unwrap();
        test_exact(0.0, 4, 0.0, variance);
        test_exact(0.3, 3, 0.63, variance);
        test_exact(1.0, 2, 0.0, variance);
    }

    #[test]
    fn test_entropy() {
        let entropy = |x: Binomial| x.entropy().unwrap();
        test_exact(0.0, 4, 0.0, entropy);
        test_absolute(0.3, 3, 1.1404671643037712668976423399228972051669206536461, 1e-15, entropy);
        test_exact(1.0, 2, 0.0, entropy);
    }

    #[test]
    fn test_skewness() {
        let skewness = |x: Binomial| x.skewness().unwrap();
        test_exact(0.0, 4, f64::INFINITY, skewness);
        test_exact(0.3, 3, 0.503952630678969636286, skewness);
        test_exact(1.0, 2, f64::NEG_INFINITY, skewness);
    }

    #[test]
    fn test_median() {
        let median = |x: Binomial| x.median();
        test_exact(0.0, 4, 0.0, median);
        test_exact(0.3, 3, 0.0, median);
        test_exact(1.0, 2, 2.0, median);
    }

    #[test]
    fn test_mode() {
        let mode = |x: Binomial| x.mode().unwrap();
        test_exact(0.0, 4, 0, mode);
        test_exact(0.3, 3, 1, mode);
        test_exact(1.0, 2, 2, mode);
    }

    #[test]
    fn test_min_max() {
        let min = |x: Binomial| x.min();
        let max = |x: Binomial| x.max();
        test_exact(0.3, 10, 0, min);
        test_exact(0.3, 10, 10, max);
    }

    #[test]
    fn test_pmf() {
        let pmf = |arg: u64| move |x: Binomial| x.pmf(arg);
        test_exact(0.0, 1, 1.0, pmf(0));
        test_exact(0.0, 1, 0.0, pmf(1));
        test_exact(0.0, 3, 1.0, pmf(0));
        test_exact(0.0, 3, 0.0, pmf(1));
        test_exact(0.0, 3, 0.0, pmf(3));
        test_exact(0.0, 10, 1.0, pmf(0));
        test_exact(0.0, 10, 0.0, pmf(1));
        test_exact(0.0, 10, 0.0, pmf(10));
        test_exact(0.3, 1, 0.69999999999999995559107901499373838305473327636719, pmf(0));
        test_exact(0.3, 1, 0.2999999999999999888977697537484345957636833190918, pmf(1));
        test_exact(0.3, 3, 0.34299999999999993471888615204079956461021032657166, pmf(0));
        test_absolute(0.3, 3, 0.44099999999999992772448109690231306411849135972008, 1e-15, pmf(1));
        test_absolute(0.3, 3, 0.026999999999999997002397833512077451789759292859569, 1e-16, pmf(3));
        test_absolute(0.3, 10, 0.02824752489999998207939855277004937778546385011091, 1e-17, pmf(0));
        test_absolute(0.3, 10, 0.12106082099999992639752977030555903089040470780077, 1e-15, pmf(1));
        test_absolute(0.3, 10, 0.0000059048999999999978147480206303047454017251032868501, 1e-20, pmf(10));
        test_exact(1.0, 1, 0.0, pmf(0));
        test_exact(1.0, 1, 1.0, pmf(1));
        test_exact(1.0, 3, 0.0, pmf(0));
        test_exact(1.0, 3, 0.0, pmf(1));
        test_exact(1.0, 3, 1.0, pmf(3));
        test_exact(1.0, 10, 0.0, pmf(0));
        test_exact(1.0, 10, 0.0, pmf(1));
        test_exact(1.0, 10, 1.0, pmf(10));
    }

    #[test]
    fn test_ln_pmf() {
        let ln_pmf = |arg: u64| move |x: Binomial| x.ln_pmf(arg);
        test_exact(0.0, 1, 0.0, ln_pmf(0));
        test_exact(0.0, 1, f64::NEG_INFINITY, ln_pmf(1));
        test_exact(0.0, 3, 0.0, ln_pmf(0));
        test_exact(0.0, 3, f64::NEG_INFINITY, ln_pmf(1));
        test_exact(0.0, 3, f64::NEG_INFINITY, ln_pmf(3));
        test_exact(0.0, 10, 0.0, ln_pmf(0));
        test_exact(0.0, 10, f64::NEG_INFINITY, ln_pmf(1));
        test_exact(0.0, 10, f64::NEG_INFINITY, ln_pmf(10));
        test_exact(0.3, 1, -0.3566749439387324423539544041072745145718090708995, ln_pmf(0));
        test_exact(0.3, 1, -1.2039728043259360296301803719337238685164245381839, ln_pmf(1));
        test_exact(0.3, 3, -1.0700248318161973270618632123218235437154272126985, ln_pmf(0));
        test_absolute(0.3, 3, -0.81871040353529122294284394322574719301255212216016, 1e-15, ln_pmf(1));
        test_absolute(0.3, 3, -3.6119184129778080888905411158011716055492736145517, 1e-15, ln_pmf(3));
        test_exact(0.3, 10, -3.566749439387324423539544041072745145718090708995, ln_pmf(0));
        test_absolute(0.3, 10, -2.1114622067804823267977785542148302920616046876506, 1e-14, ln_pmf(1));
        test_exact(0.3, 10, -12.039728043259360296301803719337238685164245381839, ln_pmf(10));
        test_exact(1.0, 1, f64::NEG_INFINITY, ln_pmf(0));
        test_exact(1.0, 1, 0.0, ln_pmf(1));
        test_exact(1.0, 3, f64::NEG_INFINITY, ln_pmf(0));
        test_exact(1.0, 3, f64::NEG_INFINITY, ln_pmf(1));
        test_exact(1.0, 3, 0.0, ln_pmf(3));
        test_exact(1.0, 10, f64::NEG_INFINITY, ln_pmf(0));
        test_exact(1.0, 10, f64::NEG_INFINITY, ln_pmf(1));
        test_exact(1.0, 10, 0.0, ln_pmf(10));
    }

    #[test]
    fn test_cdf() {
        let cdf = |arg: u64| move |x: Binomial| x.cdf(arg);
        test_exact(0.0, 1, 1.0, cdf(0));
        test_exact(0.0, 1, 1.0, cdf(1));
        test_exact(0.0, 3, 1.0, cdf(0));
        test_exact(0.0, 3, 1.0, cdf(1));
        test_exact(0.0, 3, 1.0, cdf(3));
        test_exact(0.0, 10, 1.0, cdf(0));
        test_exact(0.0, 10, 1.0, cdf(1));
        test_exact(0.0, 10, 1.0, cdf(10));
        test_absolute(0.3, 1, 0.7, 1e-15, cdf(0));
        test_exact(0.3, 1, 1.0, cdf(1));
        test_absolute(0.3, 3, 0.343, 1e-14, cdf(0));
        test_absolute(0.3, 3, 0.784, 1e-15, cdf(1));
        test_exact(0.3, 3, 1.0, cdf(3));
        test_absolute(0.3, 10, 0.0282475249, 1e-16, cdf(0));
        test_absolute(0.3, 10, 0.1493083459, 1e-14, cdf(1));
        test_exact(0.3, 10, 1.0, cdf(10));
        test_exact(1.0, 1, 0.0, cdf(0));
        test_exact(1.0, 1, 1.0, cdf(1));
        test_exact(1.0, 3, 0.0, cdf(0));
        test_exact(1.0, 3, 0.0, cdf(1));
        test_exact(1.0, 3, 1.0, cdf(3));
        test_exact(1.0, 10, 0.0, cdf(0));
        test_exact(1.0, 10, 0.0, cdf(1));
        test_exact(1.0, 10, 1.0, cdf(10));
    }

    #[test]
    fn test_sf() {
        let sf = |arg: u64| move |x: Binomial| x.sf(arg);
        test_exact(0.0, 1, 0.0, sf(0));
        test_exact(0.0, 1, 0.0, sf(1));
        test_exact(0.0, 3, 0.0, sf(0));
        test_exact(0.0, 3, 0.0, sf(1));
        test_exact(0.0, 3, 0.0, sf(3));
        test_exact(0.0, 10, 0.0, sf(0));
        test_exact(0.0, 10, 0.0, sf(1));
        test_exact(0.0, 10, 0.0, sf(10));
        test_absolute(0.3, 1, 0.3, 1e-15, sf(0));
        test_exact(0.3, 1, 0.0, sf(1));
        test_absolute(0.3, 3, 0.657, 1e-14, sf(0));
        test_absolute(0.3, 3, 0.216, 1e-15, sf(1));
        test_exact(0.3, 3, 0.0, sf(3));
        test_absolute(0.3, 10, 0.9717524751000001, 1e-16, sf(0));
        test_absolute(0.3, 10, 0.850691654100002, 1e-14, sf(1));
        test_exact(0.3, 10, 0.0, sf(10));
        test_exact(1.0, 1, 1.0, sf(0));
        test_exact(1.0, 1, 0.0, sf(1));
        test_exact(1.0, 3, 1.0, sf(0));
        test_exact(1.0, 3, 1.0, sf(1));
        test_exact(1.0, 3, 0.0, sf(3));
        test_exact(1.0, 10, 1.0, sf(0));
        test_exact(1.0, 10, 1.0, sf(1));
        test_exact(1.0, 10, 0.0, sf(10));
    }

    #[test]
    fn test_cdf_upper_bound() {
        let cdf = |arg: u64| move |x: Binomial| x.cdf(arg);
        test_exact(0.5, 3, 1.0, cdf(5));
    }

    #[test]
    fn test_sf_upper_bound() {
        let sf = |arg: u64| move |x: Binomial| x.sf(arg);
        test_exact(0.5, 3, 0.0, sf(5));
    }

    #[test]
    fn test_inverse_cdf() {
        let invcdf = |arg: f64| move |x: Binomial| x.inverse_cdf(arg);
        test_exact(0.4, 5, 2, invcdf(0.3456));

        // cases in issue #185
        test_exact(0.018, 465, 1, invcdf(3.472e-4));
        test_exact(0.5, 6, 4, invcdf(0.75));

        // case in issue #330
        test_exact(0.05, 2, 0, invcdf(0.5));      
        test_exact(0.005, 10, 0, invcdf(0.9));      
    }

    #[test]
    fn test_cdf_inverse_cdf() {
        let cdf_invcdf = |arg: u64| move |x: Binomial| x.inverse_cdf(x.cdf(arg));
        test_exact(0.3, 10, 3, cdf_invcdf(3));
        test_exact(0.3, 10, 4, cdf_invcdf(4));
        test_exact(0.5, 6, 4, cdf_invcdf(4));
    }

    #[test]
    fn test_discrete() {
        density_util::check_discrete_distribution(&create_ok(0.3, 5), 5);
        density_util::check_discrete_distribution(&create_ok(0.7, 10), 10);
    }
}
