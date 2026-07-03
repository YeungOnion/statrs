use super::Beta;
#[cfg(test)]
use super::BetaError;
use crate::distribution::{Continuous, ContinuousCDF};
use crate::statistics::*;

impl ContinuousCDF<f64, f64> for Beta {
    fn cdf(&self, x: f64) -> f64 {
        super::cdf(self.shape_a, self.shape_b, x)
    }

    fn sf(&self, x: f64) -> f64 {
        super::sf(self.shape_a, self.shape_b, x)
    }

    /// # Panics
    ///
    /// If x is not in `[0, 1]`.
    fn inverse_cdf(&self, x: f64) -> f64 {
        if !(0.0..=1.0).contains(&x) {
            panic!("x must be in [0, 1]");
        }
        super::inverse_cdf_unchecked(self.shape_a, self.shape_b, x)
    }
}

impl Min<f64> for Beta {
    fn min(&self) -> f64 {
        0.0
    }
}

impl Max<f64> for Beta {
    fn max(&self) -> f64 {
        1.0
    }
}

impl Distribution<f64> for Beta {
    fn mean(&self) -> Option<f64> {
        Some(super::mean(self.shape_a, self.shape_b))
    }

    fn variance(&self) -> Option<f64> {
        Some(super::variance(self.shape_a, self.shape_b))
    }

    fn entropy(&self) -> Option<f64> {
        Some(
            crate::function::beta::ln_beta(self.shape_a, self.shape_b)
                - (self.shape_a - 1.0) * crate::function::gamma::digamma(self.shape_a)
                - (self.shape_b - 1.0) * crate::function::gamma::digamma(self.shape_b)
                + (self.shape_a + self.shape_b - 2.0)
                    * crate::function::gamma::digamma(self.shape_a + self.shape_b),
        )
    }

    fn skewness(&self) -> Option<f64> {
        Some(
            2.0 * (self.shape_b - self.shape_a) * (self.shape_a + self.shape_b + 1.0).sqrt()
                / ((self.shape_a + self.shape_b + 2.0) * (self.shape_a * self.shape_b).sqrt()),
        )
    }
}

impl Mode<Option<f64>> for Beta {
    fn mode(&self) -> Option<f64> {
        if self.shape_a <= 1.0 || self.shape_b <= 1.0 {
            None
        } else {
            Some((self.shape_a - 1.0) / (self.shape_a + self.shape_b - 2.0))
        }
    }
}

impl Continuous<f64, f64> for Beta {
    fn pdf(&self, x: f64) -> f64 {
        super::pdf(self.shape_a, self.shape_b, x)
    }

    fn ln_pdf(&self, x: f64) -> f64 {
        super::ln_pdf(self.shape_a, self.shape_b, x)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::internal::density_util;

    crate::distribution::internal::testing_boiler!(a: f64, b: f64; Beta; BetaError);

    #[test]
    fn test_create() {
        let valid = [(1.0, 1.0), (9.0, 1.0), (5.0, 100.0)];
        for (a, b) in valid {
            create_ok(a, b);
        }
    }

    #[test]
    fn test_bad_create() {
        let invalid = [
            (0.0, 0.0),
            (0.0, 0.1),
            (1.0, 0.0),
            (0.5, f64::INFINITY),
            (f64::INFINITY, 0.5),
            (f64::NAN, 1.0),
            (1.0, f64::NAN),
            (f64::NAN, f64::NAN),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
            (f64::INFINITY, f64::INFINITY),
        ];
        for (a, b) in invalid {
            create_err(a, b);
        }
    }

    #[test]
    fn test_mean() {
        let f = |x: Beta| x.mean().unwrap();
        let test = [
            ((1.0, 1.0), 0.5),
            ((9.0, 1.0), 0.9),
            ((5.0, 100.0), 0.047619047619047619047616),
        ];
        for ((a, b), res) in test {
            test_relative(a, b, res, f);
        }
    }

    #[test]
    fn test_variance() {
        let f = |x: Beta| x.variance().unwrap();
        let test = [
            ((1.0, 1.0), 1.0 / 12.0),
            ((9.0, 1.0), 9.0 / 1100.0),
            ((5.0, 100.0), 500.0 / 1168650.0),
        ];
        for ((a, b), res) in test {
            test_relative(a, b, res, f);
        }
    }

    #[test]
    fn test_entropy() {
        let f = |x: Beta| x.entropy().unwrap();
        let test = [
            ((9.0, 1.0), -1.3083356884473304939016015),
            ((5.0, 100.0), -2.52016231876027436794592),
        ];
        for ((a, b), res) in test {
            test_relative(a, b, res, f);
        }
        test_absolute(1.0, 1.0, 0.0, 1e-14, f);
    }

    #[test]
    fn test_skewness() {
        let skewness = |x: Beta| x.skewness().unwrap();
        test_relative(1.0, 1.0, 0.0, skewness);
        test_relative(9.0, 1.0, -1.4740554623801777107177478829, skewness);
        test_relative(5.0, 100.0, 0.817594109275534303545831591, skewness);
    }

    #[test]
    fn test_mode() {
        let mode = |x: Beta| x.mode().unwrap();
        test_relative(5.0, 100.0, 0.038834951456310676243255386, mode);
    }

    #[test]
    fn test_mode_shape_a_lte_1() {
        test_none(1.0, 5.0, |dist| dist.mode());
    }

    #[test]
    fn test_mode_shape_b_lte_1() {
        test_none(5.0, 1.0, |dist| dist.mode());
    }

    #[test]
    fn test_min_max() {
        let min = |x: Beta| x.min();
        let max = |x: Beta| x.max();
        test_relative(1.0, 1.0, 0.0, min);
        test_relative(1.0, 1.0, 1.0, max);
    }

    #[test]
    fn test_pdf() {
        let x = create_ok(9.0, 1.0);
        assert_eq!(x.pdf(0.5), super::super::pdf(9.0, 1.0, 0.5));
    }

    #[test]
    fn test_ln_pdf() {
        let x = create_ok(9.0, 1.0);
        assert_eq!(x.ln_pdf(0.5), super::super::ln_pdf(9.0, 1.0, 0.5));
    }

    #[test]
    fn test_cdf() {
        let x = create_ok(9.0, 1.0);
        assert_eq!(x.cdf(0.5), super::super::cdf(9.0, 1.0, 0.5));
    }

    #[test]
    fn test_sf() {
        let x = create_ok(9.0, 1.0);
        assert_eq!(x.sf(0.5), super::super::sf(9.0, 1.0, 0.5));
    }

    #[test]
    fn test_inverse_cdf() {
        let x = create_ok(9.0, 1.0);
        assert_eq!(
            x.inverse_cdf(0.5),
            super::super::inverse_cdf_unchecked(9.0, 1.0, 0.5)
        );
    }

    #[test]
    #[should_panic]
    fn test_inverse_cdf_panics_outside_domain() {
        create_ok(9.0, 1.0).inverse_cdf(1.5);
    }

    #[test]
    fn test_continuous() {
        density_util::check_continuous_distribution(&create_ok(1.2, 3.4), 0.0, 1.0);
        density_util::check_continuous_distribution(&create_ok(4.5, 6.7), 0.0, 1.0);
    }
}
