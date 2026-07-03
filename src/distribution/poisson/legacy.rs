use super::Poisson;
#[cfg(test)]
use super::PoissonError;
use crate::distribution::{Discrete, DiscreteCDF};
use crate::statistics::*;
use core::f64;

impl DiscreteCDF<u64, f64> for Poisson {
    fn cdf(&self, x: u64) -> f64 {
        super::cdf(self.lambda, x)
    }

    fn sf(&self, x: u64) -> f64 {
        super::sf(self.lambda, x)
    }

    /// # Panics
    ///
    /// If `p` is not in `[0, 1]`.
    fn inverse_cdf(&self, p: f64) -> u64 {
        if !(0.0..=1.0).contains(&p) {
            panic!("p must be on [0, 1]");
        }
        super::inverse_cdf_unchecked(self.lambda, p)
    }
}

impl Min<u64> for Poisson {
    fn min(&self) -> u64 {
        0
    }
}

impl Max<u64> for Poisson {
    fn max(&self) -> u64 {
        u64::MAX
    }
}

impl Distribution<f64> for Poisson {
    fn mean(&self) -> Option<f64> {
        Some(super::mean(self.lambda))
    }

    fn variance(&self) -> Option<f64> {
        Some(super::variance(self.lambda))
    }

    fn entropy(&self) -> Option<f64> {
        Some(
            0.5 * (2.0 * f64::consts::PI * f64::consts::E * self.lambda).ln()
                - 1.0 / (12.0 * self.lambda)
                - 1.0 / (24.0 * self.lambda * self.lambda)
                - 19.0 / (360.0 * self.lambda * self.lambda * self.lambda),
        )
    }

    fn skewness(&self) -> Option<f64> {
        Some(1.0 / self.lambda.sqrt())
    }
}

impl Median<f64> for Poisson {
    fn median(&self) -> f64 {
        (self.lambda + 1.0 / 3.0 - 0.02 / self.lambda).floor()
    }
}

impl Mode<Option<u64>> for Poisson {
    fn mode(&self) -> Option<u64> {
        Some(self.lambda.floor() as u64)
    }
}

impl Discrete<u64, f64> for Poisson {
    fn pmf(&self, x: u64) -> f64 {
        super::pmf(self.lambda, x)
    }

    fn ln_pmf(&self, x: u64) -> f64 {
        super::ln_pmf(self.lambda, x)
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::internal::density_util;
    crate::distribution::internal::testing_boiler!(lambda: f64; Poisson; PoissonError);

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
        let x = create_ok(1.5);
        assert_eq!(x.pmf(1), super::super::pmf(1.5, 1));
    }

    #[test]
    fn test_ln_pmf() {
        let x = create_ok(1.5);
        assert_eq!(x.ln_pmf(1), super::super::ln_pmf(1.5, 1));
    }

    #[test]
    fn test_cdf() {
        let x = create_ok(1.5);
        assert_eq!(x.cdf(1), super::super::cdf(1.5, 1));
    }

    #[test]
    fn test_sf() {
        let x = create_ok(1.5);
        assert_eq!(x.sf(1), super::super::sf(1.5, 1));
    }

    #[test]
    fn test_discrete() {
        density_util::check_discrete_distribution(&create_ok(0.3), 10);
        density_util::check_discrete_distribution(&create_ok(4.5), 30);
    }

    #[test]
    fn test_inverse_cdf() {
        let x = create_ok(1.5);
        assert_eq!(x.inverse_cdf(0.5), super::super::inverse_cdf_unchecked(1.5, 0.5));
    }

    #[test]
    #[should_panic]
    fn test_inverse_cdf_panics_outside_domain() {
        create_ok(1.5).inverse_cdf(1.5);
    }
}
