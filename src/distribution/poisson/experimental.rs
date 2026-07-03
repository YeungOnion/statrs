use super::Poisson;
use crate::experimental_api::{
    ClosedFormCdf, InvalidVariate, InverseCdf, InverseCdfError, Mean, Pmf, Probability,
    ProbabilityMass, TryVariate, Variance, Variate,
};

impl crate::Sealed for Poisson {}

impl TryVariate for Poisson {
    type Repr = u64;

    fn try_variate(&self, x: u64) -> Result<Variate<Self, u64>, InvalidVariate<u64>> {
        Ok(Variate::new(x))
    }
}

impl ClosedFormCdf for Poisson {
    fn cdf(&self, x: Variate<Self, u64>) -> Probability {
        Probability::new(super::cdf(self.lambda, x.into_inner()))
            .expect("Poisson CDF is always in [0, 1]")
    }

    fn sf(&self, x: Variate<Self, u64>) -> Probability {
        let sf = crate::function::gamma::gamma_lr(x.into_inner() as f64 + 1.0, self.lambda);
        Probability::new(sf).expect("Poisson SF is always in [0, 1]")
    }
}

impl InverseCdf for Poisson {
    fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, u64>, InverseCdfError> {
        let x = super::inverse_cdf_unchecked(self.lambda, p.into_inner());
        self.try_variate(x)
            .map_err(|_| InverseCdfError::OutOfSupport)
    }
}

impl Pmf for Poisson {
    fn pmf(&self, x: Variate<Self, u64>) -> ProbabilityMass {
        ProbabilityMass::new(super::pmf(self.lambda, x.into_inner()))
            .expect("Poisson PMF is always in [0, 1]")
    }
}

impl Mean for Poisson {
    type Output = f64;
    fn mean(&self) -> f64 {
        super::mean(self.lambda)
    }
}

impl Variance for Poisson {
    type Output = f64;
    fn variance(&self) -> f64 {
        super::variance(self.lambda)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn try_variate_always_succeeds() {
        let d = Poisson::new(1.5).unwrap();
        assert!(d.try_variate(0).is_ok());
        assert!(d.try_variate(u64::MAX).is_ok());
    }

    #[test]
    fn pmf_delegates_to_shared_fn() {
        let d = Poisson::new(5.4).unwrap();
        let x = d.try_variate(10).unwrap();
        assert_eq!(d.pmf(x).into_inner(), super::super::pmf(5.4, 10));
    }

    #[test]
    fn cdf_delegates_to_shared_fn() {
        let d = Poisson::new(5.4).unwrap();
        let x = d.try_variate(10).unwrap();
        assert_eq!(d.cdf(x).into_inner(), super::super::cdf(5.4, 10));
    }

    #[test]
    fn sf_complements_cdf() {
        let d = Poisson::new(5.4).unwrap();
        let x = d.try_variate(10).unwrap();
        crate::prec::assert_abs_diff_eq!(
            d.sf(x).into_inner(),
            1.0 - d.cdf(x).into_inner(),
            epsilon = 1e-10
        );
    }

    #[test]
    fn mean_matches_shared_fn() {
        let d = Poisson::new(5.4).unwrap();
        assert_eq!(d.mean(), super::super::mean(5.4));
    }

    #[test]
    fn variance_matches_shared_fn() {
        let d = Poisson::new(5.4).unwrap();
        assert_eq!(d.variance(), super::super::variance(5.4));
    }

    #[test]
    fn inverse_cdf_roundtrips_through_cdf() {
        let d = Poisson::new(1.5).unwrap();
        let p = Probability::new(0.5).unwrap();
        let x = d.inverse_cdf(p).unwrap();
        assert!(d.cdf(x).into_inner() >= 0.5);
    }
}
