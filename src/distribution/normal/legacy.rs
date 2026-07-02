use super::{Normal, NormalError, cdf_unchecked, ln_pdf_unchecked, pdf_unchecked, sf_unchecked};
use crate::consts;
use crate::distribution::{Continuous, ContinuousCDF};
use crate::statistics::*;
use core::f64;

impl ContinuousCDF<f64, f64> for Normal {
    fn cdf(&self, x: f64) -> f64 {
        cdf_unchecked(x, self.mean, self.std_dev)
    }

    fn sf(&self, x: f64) -> f64 {
        sf_unchecked(x, self.mean, self.std_dev)
    }

    /// # Panics
    ///
    /// If `x < 0.0` or `x > 1.0`
    fn inverse_cdf(&self, x: f64) -> f64 {
        if !(0.0..=1.0).contains(&x) {
            panic!("x must be in [0, 1]");
        }
        super::inverse_cdf_unchecked(x, self.mean, self.std_dev)
    }
}

impl Min<f64> for Normal {
    fn min(&self) -> f64 {
        f64::NEG_INFINITY
    }
}

impl Max<f64> for Normal {
    fn max(&self) -> f64 {
        f64::INFINITY
    }
}

impl Distribution<f64> for Normal {
    fn mean(&self) -> Option<f64> {
        Some(self.mean)
    }

    fn variance(&self) -> Option<f64> {
        Some(self.std_dev * self.std_dev)
    }

    fn std_dev(&self) -> Option<f64> {
        Some(self.std_dev)
    }

    fn entropy(&self) -> Option<f64> {
        Some(self.std_dev.ln() + consts::LN_SQRT_2PIE)
    }

    fn skewness(&self) -> Option<f64> {
        Some(0.0)
    }
}

impl Median<f64> for Normal {
    fn median(&self) -> f64 {
        self.mean
    }
}

impl Mode<Option<f64>> for Normal {
    fn mode(&self) -> Option<f64> {
        Some(self.mean)
    }
}

impl Continuous<f64, f64> for Normal {
    fn pdf(&self, x: f64) -> f64 {
        pdf_unchecked(x, self.mean, self.std_dev)
    }

    fn ln_pdf(&self, x: f64) -> f64 {
        ln_pdf_unchecked(x, self.mean, self.std_dev)
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::prec;
    use crate::distribution::internal::density_util;
    crate::distribution::internal::testing_boiler!(mean: f64, std_dev: f64; Normal; NormalError);

    #[test]
    fn test_create() {
        create_ok(10.0, 0.1);
        create_ok(-5.0, 1.0);
        create_ok(0.0, 10.0);
        create_ok(10.0, 100.0);
        create_ok(-5.0, f64::INFINITY);
    }

    #[test]
    fn test_bad_create() {
        test_create_err(f64::NAN, 1.0, NormalError::MeanInvalid);
        test_create_err(1.0, f64::NAN, NormalError::StandardDeviationInvalid);
        create_err(0.0, 0.0);
        create_err(f64::NAN, f64::NAN);
        create_err(1.0, -1.0);
    }

    #[test]
    fn test_variance() {
        let variance = |x: Normal| x.variance().unwrap();
        test_exact(0.0, 0.1, 0.1 * 0.1, variance);
        test_exact(0.0, 1.0, 1.0, variance);
        test_exact(0.0, 10.0, 100.0, variance);
        test_exact(0.0, f64::INFINITY, f64::INFINITY, variance);
    }

    #[test]
    fn test_entropy() {
        let entropy = |x: Normal| x.entropy().unwrap();
        test_absolute(0.0, 0.1, -0.8836465597893729422377, 1e-15, entropy);
        test_exact(0.0, 1.0, 1.41893853320467274178, entropy);
        test_exact(0.0, 10.0, 3.721523626198718425798, entropy);
        test_exact(0.0, f64::INFINITY, f64::INFINITY, entropy);
    }

    #[test]
    fn test_skewness() {
        let skewness = |x: Normal| x.skewness().unwrap();
        test_exact(0.0, 0.1, 0.0, skewness);
        test_exact(4.0, 1.0, 0.0, skewness);
        test_exact(0.3, 10.0, 0.0, skewness);
        test_exact(0.0, f64::INFINITY, 0.0, skewness);
    }

    #[test]
    fn test_mode() {
        let mode = |x: Normal| x.mode().unwrap();
        test_exact(-0.0, 1.0, 0.0, mode);
        test_exact(0.0, 1.0, 0.0, mode);
        test_exact(0.1, 1.0, 0.1, mode);
        test_exact(1.0, 1.0, 1.0, mode);
        test_exact(-10.0, 1.0, -10.0, mode);
        test_exact(f64::INFINITY, 1.0, f64::INFINITY, mode);
    }

    #[test]
    fn test_median() {
        let median = |x: Normal| x.median();
        test_exact(-0.0, 1.0, 0.0, median);
        test_exact(0.0, 1.0, 0.0, median);
        test_exact(0.1, 1.0, 0.1, median);
        test_exact(1.0, 1.0, 1.0, median);
        test_exact(-0.0, 1.0, -0.0, median);
        test_exact(f64::INFINITY, 1.0, f64::INFINITY, median);
    }

    #[test]
    fn test_min_max() {
        let min = |x: Normal| x.min();
        let max = |x: Normal| x.max();
        test_exact(0.0, 0.1, f64::NEG_INFINITY, min);
        test_exact(-3.0, 10.0, f64::NEG_INFINITY, min);
        test_exact(0.0, 0.1, f64::INFINITY, max);
        test_exact(-3.0, 10.0, f64::INFINITY, max);
    }

    #[test]
    fn test_pdf() {
        let x = create_ok(10.0, 0.1);
        assert_eq!(x.pdf(10.0), pdf_unchecked(10.0, 10.0, 0.1));
    }

    #[test]
    fn test_ln_pdf() {
        let x = create_ok(10.0, 0.1);
        assert_eq!(x.ln_pdf(10.0), ln_pdf_unchecked(10.0, 10.0, 0.1));
    }

    #[test]
    fn test_cdf() {
        let x = create_ok(5.0, 2.0);
        assert_eq!(x.cdf(4.0), cdf_unchecked(4.0, 5.0, 2.0));
    }

    #[test]
    fn test_sf() {
        let x = create_ok(5.0, 2.0);
        assert_eq!(x.sf(4.0), sf_unchecked(4.0, 5.0, 2.0));
    }

    #[test]
    fn test_continuous() {
        density_util::check_continuous_distribution(&create_ok(0.0, 1.0), -10.0, 10.0);
        density_util::check_continuous_distribution(&create_ok(20.0, 0.5), 10.0, 30.0);
    }

    #[test]
    fn test_inverse_cdf() {
        let x = create_ok(5.0, 2.0);
        assert_eq!(x.inverse_cdf(0.5), super::super::inverse_cdf_unchecked(0.5, 5.0, 2.0));
    }

    #[test]
    #[should_panic]
    fn test_inverse_cdf_panics_outside_domain() {
        create_ok(5.0, 2.0).inverse_cdf(1.5);
    }

    #[test]
    fn test_default() {
        let n = Normal::default();
        let n_mean = n.mean().unwrap();
        let n_std = n.std_dev().unwrap();
        prec::assert_abs_diff_eq!(n_mean, 0.0, epsilon = 1e-15);
        prec::assert_abs_diff_eq!(n_std, 1.0, epsilon = 1e-15);
    }
}
