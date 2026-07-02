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
#[cfg(feature = "experimental_api")]
mod experimental;

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

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: f64, expect: f64, epsilon: f64) {
        if expect.is_infinite() && actual == expect {
            return;
        }
        crate::prec::assert_abs_diff_eq!(actual, expect, epsilon = epsilon);
    }

    #[test]
    fn pdf_unchecked_matches_reference_values() {
        let test = [
            (10.0, 0.1, 8.5, 5.530709549844416159162E-49, 1e-64),
            (10.0, 0.1, 9.8, 0.5399096651318805195056, 1e-14),
            (10.0, 0.1, 10.0, 3.989422804014326779399, 1e-15),
            (10.0, 0.1, 10.2, 0.5399096651318805195056, 1e-14),
            (10.0, 0.1, 11.5, 5.530709549844416159162E-49, 1e-64),
            (-5.0, 1.0, -10.0, 1.486719514734297707908E-6, 0.0),
            (-5.0, 1.0, -7.5, 0.01752830049356853736216, 0.0),
            (-5.0, 1.0, -5.0, 0.3989422804014326779399, 1e-16),
            (-5.0, 1.0, -2.5, 0.01752830049356853736216, 0.0),
            (-5.0, 1.0, 0.0, 1.486719514734297707908E-6, 0.0),
            (0.0, 10.0, -5.0, 0.03520653267642994777747, 0.0),
            (0.0, 10.0, -2.5, 0.03866681168028492069412, 1e-17),
            (0.0, 10.0, 0.0, 0.03989422804014326779399, 1e-17),
            (0.0, 10.0, 2.5, 0.03866681168028492069412, 1e-17),
            (0.0, 10.0, 5.0, 0.03520653267642994777747, 0.0),
            (10.0, 100.0, -200.0, 4.398359598042719404845E-4, 1e-19),
            (10.0, 100.0, -100.0, 0.002178521770325505313831, 0.0),
            (10.0, 100.0, 0.0, 0.003969525474770117655105, 0.0),
            (10.0, 100.0, 100.0, 0.002660852498987548218204, 1e-18),
            (10.0, 100.0, 200.0, 6.561581477467659126534E-4, 0.0),
            (-5.0, f64::INFINITY, -5.0, 0.0, 0.0),
            (-5.0, f64::INFINITY, 0.0, 0.0, 0.0),
            (-5.0, f64::INFINITY, 100.0, 0.0, 0.0),
        ];
        for (mean, std_dev, x, expect, eps) in test {
            assert_close(pdf_unchecked(x, mean, std_dev), expect, eps);
        }
    }

    #[test]
    fn ln_pdf_unchecked_matches_reference_values() {
        let test = [
            (10.0, 0.1, 8.5, (5.530709549844416159162E-49f64).ln(), 1e-13),
            (10.0, 0.1, 9.8, (0.5399096651318805195056f64).ln(), 1e-13),
            (10.0, 0.1, 10.0, (3.989422804014326779399f64).ln(), 1e-15),
            (10.0, 0.1, 10.2, (0.5399096651318805195056f64).ln(), 1e-13),
            (10.0, 0.1, 11.5, (5.530709549844416159162E-49f64).ln(), 1e-13),
            (-5.0, 1.0, -10.0, (1.486719514734297707908E-6f64).ln(), 0.0),
            (-5.0, 1.0, -7.5, (0.01752830049356853736216f64).ln(), 0.0),
            (-5.0, 1.0, -5.0, (0.3989422804014326779399f64).ln(), 1e-15),
            (-5.0, 1.0, -2.5, (0.01752830049356853736216f64).ln(), 0.0),
            (-5.0, 1.0, 0.0, (1.486719514734297707908E-6f64).ln(), 0.0),
            (0.0, 10.0, -5.0, (0.03520653267642994777747f64).ln(), 0.0),
            (0.0, 10.0, -2.5, (0.03866681168028492069412f64).ln(), 0.0),
            (0.0, 10.0, 0.0, (0.03989422804014326779399f64).ln(), 0.0),
            (0.0, 10.0, 2.5, (0.03866681168028492069412f64).ln(), 0.0),
            (0.0, 10.0, 5.0, (0.03520653267642994777747f64).ln(), 0.0),
            (10.0, 100.0, -200.0, (4.398359598042719404845E-4f64).ln(), 0.0),
            (10.0, 100.0, -100.0, (0.002178521770325505313831f64).ln(), 0.0),
            (10.0, 100.0, 0.0, (0.003969525474770117655105f64).ln(), 1e-15),
            (10.0, 100.0, 100.0, (0.002660852498987548218204f64).ln(), 1e-15),
            (10.0, 100.0, 200.0, (6.561581477467659126534E-4f64).ln(), 1e-15),
            (-5.0, f64::INFINITY, -5.0, f64::NEG_INFINITY, 0.0),
            (-5.0, f64::INFINITY, 0.0, f64::NEG_INFINITY, 0.0),
            (-5.0, f64::INFINITY, 100.0, f64::NEG_INFINITY, 0.0),
        ];
        for (mean, std_dev, x, expect, eps) in test {
            assert_close(ln_pdf_unchecked(x, mean, std_dev), expect, eps);
        }
    }

    #[test]
    fn cdf_unchecked_matches_reference_values() {
        let test = [
            (5.0, 2.0, f64::NEG_INFINITY, 0.0, 0.0),
            (5.0, 2.0, -5.0, 0.0000002866515718, 1e-16),
            (5.0, 2.0, -2.0, 0.0002326290790, 1e-13),
            (5.0, 2.0, 0.0, 0.006209665325, 1e-12),
            (5.0, 2.0, 4.0, 0.30853753872598689636229538939166226011639782444542207, 0.0),
            (5.0, 2.0, 5.0, 0.5, 0.0),
            (5.0, 2.0, 6.0, 0.69146246127401310363770461060833773988360217555457859, 0.0),
            (5.0, 2.0, 10.0, 0.993790334674, 1e-12),
        ];
        for (mean, std_dev, x, expect, eps) in test {
            assert_close(cdf_unchecked(x, mean, std_dev), expect, eps);
        }
    }

    #[test]
    fn sf_unchecked_matches_reference_values() {
        let test = [
            (5.0, 2.0, f64::NEG_INFINITY, 1.0, 0.0),
            (5.0, 2.0, -5.0, 0.9999997133484281, 1e-16),
            (5.0, 2.0, -2.0, 0.9997673709209455, 1e-13),
            (5.0, 2.0, 0.0, 0.9937903346744879, 1e-12),
            (5.0, 2.0, 4.0, 0.6914624612740131, 0.0),
            (5.0, 2.0, 5.0, 0.5, 0.0),
            (5.0, 2.0, 6.0, 0.3085375387259869, 0.0),
            (5.0, 2.0, 10.0, 0.006209665325512148, 1e-12),
        ];
        for (mean, std_dev, x, expect, eps) in test {
            assert_close(sf_unchecked(x, mean, std_dev), expect, eps);
        }
    }

    #[test]
    fn inverse_cdf_unchecked_matches_reference_values() {
        let test = [
            (5.0, 2.0, 0.0, f64::NEG_INFINITY, 0.0),
            (5.0, 2.0, 0.00000028665157187919391167375233287464535385442301361187883, -5.0, 1e-14),
            (5.0, 2.0, 0.0002326290790355250363499258867279847735487493358890356, -2.0, 1e-14),
            (5.0, 2.0, 0.0062096653257761351669781045741922211278977469230927036, -0.0, 1e-14),
            (5.0, 2.0, 0.30853753872598689636229538939166226011639782444542207, 4.0, 1e-14),
            (5.0, 2.0, 0.5, 5.0, 1e-14),
            (5.0, 2.0, 0.69146246127401310363770461060833773988360217555457859, 6.0, 1e-14),
            (5.0, 2.0, 0.9937903346742238648330218954258077788721022530769078, 10.0, 1e-14),
            (5.0, 2.0, 1.0, f64::INFINITY, 0.0),
        ];
        for (mean, std_dev, p, expect, eps) in test {
            assert_close(inverse_cdf_unchecked(p, mean, std_dev), expect, eps);
        }
    }

    #[test]
    fn inverse_cdf_unchecked_roundtrips_through_cdf_unchecked() {
        let test = [(5.0, 2.0, 0.3), (0.0, 1.0, 0.975)];
        for (mean, std_dev, p) in test {
            let x = inverse_cdf_unchecked(p, mean, std_dev);
            crate::prec::assert_abs_diff_eq!(cdf_unchecked(x, mean, std_dev), p, epsilon = 1e-10);
        }
    }
}

