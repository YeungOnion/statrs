use super::Beta;
use crate::experimental_api::{
    ClosedFormCdf, InvalidVariate, InverseCdf, InverseCdfError, Mean, Pdf, Probability,
    ProbabilityDensity, TryVariate, Variance, Variate,
};

impl crate::Sealed for Beta {}

impl TryVariate for Beta {
    type Repr = f64;

    fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
        if x.is_finite() && (0.0..=1.0).contains(&x) {
            Ok(Variate::new(x))
        } else {
            Err(InvalidVariate(x))
        }
    }
}

impl ClosedFormCdf for Beta {
    fn cdf(&self, x: Variate<Self, f64>) -> Probability {
        Probability::new(super::cdf(self.shape_a, self.shape_b, x.into_inner()))
            .expect("Beta CDF is always in [0, 1]")
    }
}

impl InverseCdf for Beta {
    fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, f64>, InverseCdfError> {
        let x = super::inverse_cdf_unchecked(self.shape_a, self.shape_b, p.into_inner());
        self.try_variate(x)
            .map_err(|_| InverseCdfError::OutOfSupport)
    }
}

impl Pdf for Beta {
    fn pdf(&self, x: Variate<Self, f64>) -> ProbabilityDensity {
        ProbabilityDensity::new(super::pdf(self.shape_a, self.shape_b, x.into_inner()))
            .expect("Beta PDF is always non-negative")
    }
}

impl Mean for Beta {
    type Output = f64;
    fn mean(&self) -> f64 {
        super::mean(self.shape_a, self.shape_b)
    }
}

impl Variance for Beta {
    type Output = f64;
    fn variance(&self) -> f64 {
        super::variance(self.shape_a, self.shape_b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::distribution::BetaError;

    #[test]
    fn try_variate_accepts_interior() {
        assert!(Beta::new(2.0, 2.0).unwrap().try_variate(0.5).is_ok());
    }

    #[test]
    fn try_variate_accepts_bounds() {
        let d = Beta::new(2.0, 2.0).unwrap();
        assert!(d.try_variate(0.0).is_ok());
        assert!(d.try_variate(1.0).is_ok());
    }

    #[test]
    fn try_variate_rejects_out_of_range() {
        let d = Beta::new(2.0, 2.0).unwrap();
        assert!(d.try_variate(-0.001).is_err());
        assert!(d.try_variate(1.001).is_err());
    }

    #[test]
    fn try_variate_rejects_nan() {
        assert!(Beta::new(2.0, 2.0).unwrap().try_variate(f64::NAN).is_err());
    }

    #[test]
    fn cdf_delegates_to_shared_fn() {
        let d = Beta::new(9.0, 1.0).unwrap();
        let x = d.try_variate(0.5).unwrap();
        assert_eq!(d.cdf(x).into_inner(), super::super::cdf(9.0, 1.0, 0.5));
    }

    #[test]
    fn pdf_delegates_to_shared_fn() {
        let d = Beta::new(9.0, 1.0).unwrap();
        let x = d.try_variate(0.5).unwrap();
        assert_eq!(d.pdf(x).into_inner(), super::super::pdf(9.0, 1.0, 0.5));
    }

    #[test]
    fn mean_matches_shared_fn() {
        let d = Beta::new(5.0, 100.0).unwrap();
        assert_eq!(d.mean(), super::super::mean(5.0, 100.0));
    }

    #[test]
    fn variance_matches_shared_fn() {
        let d = Beta::new(5.0, 100.0).unwrap();
        assert_eq!(d.variance(), super::super::variance(5.0, 100.0));
    }

    #[test]
    fn inverse_cdf_roundtrips_through_cdf() {
        let d = Beta::new(9.0, 1.0).unwrap();
        let p = Probability::new(0.5).unwrap();
        let x = d.inverse_cdf(p).unwrap();
        let back = d.cdf(x);
        crate::prec::assert_abs_diff_eq!(back.into_inner(), 0.5, epsilon = 1e-10);
    }

    #[test]
    fn inverse_cdf_rejects_p_outside_domain_at_construction() {
        assert!(Probability::new(1.5).is_err());
    }

    #[test]
    fn error_type_still_send_sync() {
        fn assert_sync_send<T: Sync + Send>() {}
        assert_sync_send::<BetaError>();
    }
}
