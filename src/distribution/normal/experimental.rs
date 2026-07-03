use super::{Normal, cdf_unchecked, inverse_cdf_unchecked, pdf_unchecked};
use crate::experimental_api::{
    ClosedFormCdf, InvalidVariate, InverseCdf, InverseCdfError, Mean, Pdf, Probability,
    ProbabilityDensity, TryVariate, Variance, Variate,
};

impl crate::Sealed for Normal {}

impl TryVariate for Normal {
    type Repr = f64;

    fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
        if x.is_finite() {
            Ok(Variate::new(x))
        } else {
            Err(InvalidVariate(x))
        }
    }
}

impl ClosedFormCdf for Normal {
    fn cdf(&self, x: Variate<Self, f64>) -> Probability {
        Probability::new(cdf_unchecked(x.into_inner(), self.mean, self.std_dev))
            .expect("Normal CDF is always in [0, 1]")
    }
}

impl InverseCdf for Normal {
    fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, f64>, InverseCdfError> {
        let x = inverse_cdf_unchecked(p.into_inner(), self.mean, self.std_dev);
        self.try_variate(x)
            .map_err(|_| InverseCdfError::OutOfSupport)
    }
}

impl Pdf for Normal {
    fn pdf(&self, x: Variate<Self, f64>) -> ProbabilityDensity {
        ProbabilityDensity::new(pdf_unchecked(x.into_inner(), self.mean, self.std_dev))
            .expect("Normal PDF is always non-negative")
    }
}

impl Mean for Normal {
    type Output = f64;
    fn mean(&self) -> f64 {
        self.mean
    }
}

impl Variance for Normal {
    type Output = f64;
    fn variance(&self) -> f64 {
        self.std_dev * self.std_dev
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_variate_accepts_any_finite_value() {
        let d = Normal::new(0.0, 1.0).unwrap();
        assert!(d.try_variate(-1e300).is_ok());
        assert!(d.try_variate(1e300).is_ok());
    }

    #[test]
    fn try_variate_rejects_nan_and_infinite() {
        let d = Normal::new(0.0, 1.0).unwrap();
        assert!(d.try_variate(f64::NAN).is_err());
        assert!(d.try_variate(f64::INFINITY).is_err());
        assert!(d.try_variate(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn pdf_delegates_to_shared_fn() {
        let d = Normal::new(0.0, 1.0).unwrap();
        let x = d.try_variate(1.0).unwrap();
        assert_eq!(d.pdf(x).into_inner(), pdf_unchecked(1.0, 0.0, 1.0));
    }

    #[test]
    fn cdf_delegates_to_shared_fn() {
        let d = Normal::new(5.0, 2.0).unwrap();
        let x = d.try_variate(6.0).unwrap();
        assert_eq!(d.cdf(x).into_inner(), cdf_unchecked(6.0, 5.0, 2.0));
    }

    #[test]
    fn mean_and_variance_read_fields_directly() {
        let d = Normal::new(5.0, 2.0).unwrap();
        assert_eq!(d.mean(), 5.0);
        assert_eq!(d.variance(), 4.0);
    }

    #[test]
    fn inverse_cdf_roundtrips_through_cdf() {
        let d = Normal::new(5.0, 2.0).unwrap();
        let p = Probability::new(0.5).unwrap();
        let x = d.inverse_cdf(p).unwrap();
        crate::prec::assert_abs_diff_eq!(d.cdf(x).into_inner(), 0.5, epsilon = 1e-10);
    }
}
