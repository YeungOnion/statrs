use crate::function::{factorial, gamma};
use core::f64;

/// Implements the [Poisson](https://en.wikipedia.org/wiki/Poisson_distribution)
/// distribution
///
/// # Examples
///
/// ```
/// use statrs::distribution::Poisson;
///
/// let n = Poisson::new(1.0).unwrap();
///
/// #[cfg(feature = "experimental_api")]
/// {
///     use approx::assert_abs_diff_eq;
///     use statrs::experimental_api::{Mean, Pmf, TryVariate};
///
///     assert_eq!(n.mean(), 1.0);
///     let x = n.try_variate(1).unwrap();
///     assert_abs_diff_eq!(n.pmf(x).into_inner(), 0.367879441171442, epsilon = 1e-15);
/// }
///
/// #[cfg(not(feature = "experimental_api"))]
/// {
///     use approx::assert_abs_diff_eq;
///     use statrs::distribution::Discrete;
///     use statrs::statistics::Distribution;
///
///     assert_eq!(n.mean().unwrap(), 1.0);
///     assert_abs_diff_eq!(n.pmf(1), 0.367879441171442, epsilon = 1e-15);
/// }
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

#[cfg(not(feature = "experimental_api"))]
mod legacy;

pub(crate) fn pmf(lambda: f64, x: u64) -> f64 {
    (-lambda + x as f64 * lambda.ln() - factorial::ln_factorial(x)).exp()
}

pub(crate) fn ln_pmf(lambda: f64, x: u64) -> f64 {
    -lambda + x as f64 * lambda.ln() - factorial::ln_factorial(x)
}

pub(crate) fn cdf(lambda: f64, x: u64) -> f64 {
    gamma::gamma_ur(x as f64 + 1.0, lambda)
}

pub(crate) fn sf(lambda: f64, x: u64) -> f64 {
    gamma::gamma_lr(x as f64 + 1.0, lambda)
}

pub(crate) fn inverse_cdf_unchecked(lambda: f64, p: f64) -> u64 {
    if p <= cdf(lambda, 0) {
        return 0;
    } else if p == 1.0 {
        return u64::MAX;
    }
    let mut ub: u64 = 2;
    while cdf(lambda, ub) < p {
        ub *= 2;
    }
    crate::distribution::internal::integral_bisection_search(|x: &u64| cdf(lambda, *x), p, 0u64, ub)
        .unwrap()
}

pub(crate) fn mean(lambda: f64) -> f64 {
    lambda
}

pub(crate) fn variance(lambda: f64) -> f64 {
    lambda
}

impl core::fmt::Display for Poisson {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Pois({})", self.lambda)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl ::rand::distr::Distribution<u64> for Poisson {
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
impl ::rand::distr::Distribution<f64> for Poisson {
    /// Generates one sample from the Poisson distribution either by
    /// Knuth's method if lambda < 30.0 or Rejection method PA by
    /// A. C. Atkinson from the Journal of the Royal Statistical Society
    /// Series C (Applied Statistics) Vol. 28 No. 1. (1979) pp. 29 - 35
    /// otherwise
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        sample_unchecked(rng, self.lambda)
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
        let mut product: f64 = rng.random();
        while product >= limit {
            count += 1.0;
            product *= rng.random::<f64>();
        }
        count
    } else {
        let c = 0.767 - 3.36 / lambda;
        let beta = f64::consts::PI / (3.0 * lambda).sqrt();
        let alpha = beta * lambda;
        let k = c.ln() - lambda - beta.ln();

        loop {
            let u: f64 = rng.random();
            let x = (alpha - ((1.0 - u) / u).ln()) / beta;
            let n = (x + 0.5).floor();
            if n < 0.0 {
                continue;
            }

            let v: f64 = rng.random();
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

