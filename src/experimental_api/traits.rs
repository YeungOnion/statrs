use crate::experimental_api::bisect::{BisectDomain, PartitionSpace};
use crate::experimental_api::types::{
    CdfError, InverseCdfError, Probability, ProbabilityDensity, ProbabilityMass,
};

pub use crate::experimental_api::bisect::DEFAULT_MAX_ITER;

pub trait Cdf<K> {
    type Space: PartitionSpace<Point = K>;

    fn cdf(&self, x: K) -> Result<Probability, CdfError>;
    fn cdf_domain(&self) -> Self::Space;

    fn sf(&self, x: K) -> Result<Probability, CdfError>
    where
        K: Copy,
    {
        let p = self.cdf(x)?.into_inner();
        Probability::new(1.0 - p).map_err(|_| CdfError::InvalidInput)
    }

    /// Inverse CDF via bisection. Override with a closed-form implementation where available.
    fn inverse_cdf(&self, p: Probability) -> Result<K, InverseCdfError>
    where
        Self: Sized,
        K: BisectDomain,
        Self::Space: PartitionSpace<Point = K, Cut = K>,
    {
        K::run_bisect(self, p, self.cdf_domain()).ok_or(InverseCdfError::NoConvergence)
    }
}

pub trait Pdf<K> {
    fn pdf(&self, x: K) -> Result<ProbabilityDensity, CdfError>;
    fn ln_pdf(&self, x: K) -> Result<f64, CdfError>;
}

pub trait Pmf<K> {
    fn pmf(&self, x: K) -> Result<ProbabilityMass, CdfError>;
    fn ln_pmf(&self, x: K) -> Result<f64, CdfError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_api::bisect::Interval;
    use crate::experimental_api::types::{
        CdfError, Probability, ProbabilityDensity, ProbabilityMass,
    };

    // Uniform distribution on [0.0, 1.0]: cdf(x) = clamp(x, 0, 1)
    struct UnitUniform;

    impl Cdf<f64> for UnitUniform {
        type Space = Interval<f64>;

        fn cdf(&self, x: f64) -> Result<Probability, CdfError> {
            let p = x.clamp(0.0, 1.0);
            Probability::new(p).map_err(|_| CdfError::InvalidInput)
        }

        fn cdf_domain(&self) -> Interval<f64> {
            Interval { lo: 0.0, hi: 1.0 }
        }
    }

    impl Pdf<f64> for UnitUniform {
        fn pdf(&self, x: f64) -> Result<ProbabilityDensity, CdfError> {
            let d = if (0.0..=1.0).contains(&x) { 1.0 } else { 0.0 };
            ProbabilityDensity::new(d).map_err(|_| CdfError::InvalidInput)
        }

        fn ln_pdf(&self, x: f64) -> Result<f64, CdfError> {
            Ok(self.pdf(x)?.into_inner().ln())
        }
    }

    // Discrete uniform on {0, 1, ..., 9}: cdf(k) = (k + 1) / 10
    struct TenPoint;

    impl Cdf<u64> for TenPoint {
        type Space = Interval<u64>;

        fn cdf(&self, x: u64) -> Result<Probability, CdfError> {
            let p = (x + 1).min(10) as f64 / 10.0;
            Probability::new(p).map_err(|_| CdfError::InvalidInput)
        }

        fn cdf_domain(&self) -> Interval<u64> {
            Interval { lo: 0, hi: 10 }
        }
    }

    impl Pmf<u64> for TenPoint {
        fn pmf(&self, _x: u64) -> Result<ProbabilityMass, CdfError> {
            ProbabilityMass::new(0.1).map_err(|_| CdfError::InvalidInput)
        }

        fn ln_pmf(&self, x: u64) -> Result<f64, CdfError> {
            Ok(self.pmf(x)?.into_inner().ln())
        }
    }

    #[test]
    fn cdf_evaluates_correctly() {
        let d = UnitUniform;
        let p = d.cdf(0.5).unwrap();
        assert!((p.into_inner() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn sf_is_complement_of_cdf() {
        let d = UnitUniform;
        let cdf = d.cdf(0.3).unwrap().into_inner();
        let sf = d.sf(0.3).unwrap().into_inner();
        assert!((cdf + sf - 1.0).abs() < 1e-12);
    }

    #[test]
    fn inverse_cdf_roundtrips_continuous() {
        let d = UnitUniform;
        for v in [0.1, 0.25, 0.5, 0.75, 0.9] {
            let p = Probability::new(v).unwrap();
            let x = d.inverse_cdf(p).unwrap();
            assert!(
                (x - v).abs() < 1e-6,
                "roundtrip failed for p={v}: got x={x}"
            );
        }
    }

    #[test]
    fn inverse_cdf_roundtrips_discrete() {
        let d = TenPoint;
        // cdf(k) = (k+1)/10, so inverse_cdf(0.35) should give k=3 (cdf(3)=0.4 >= 0.35)
        let p = Probability::new(0.35).unwrap();
        let k = d.inverse_cdf(p).unwrap();
        assert_eq!(k, 3);
    }

    #[test]
    fn pdf_returns_density() {
        let d = UnitUniform;
        let density = d.pdf(0.5).unwrap().into_inner();
        assert!((density - 1.0).abs() < 1e-12);
    }

    #[test]
    fn pmf_returns_mass() {
        let d = TenPoint;
        let mass = d.pmf(3).unwrap().into_inner();
        assert!((mass - 0.1).abs() < 1e-12);
    }

    #[test]
    fn pdf_outside_support_is_zero() {
        let d = UnitUniform;
        assert_eq!(d.pdf(-1.0).unwrap().into_inner(), 0.0);
        assert_eq!(d.pdf(2.0).unwrap().into_inner(), 0.0);
    }

    #[test]
    fn inverse_cdf_boundary_probabilities_continuous() {
        let d = UnitUniform;
        let p0 = Probability::new(0.0).unwrap();
        let p1 = Probability::new(1.0).unwrap();
        let x0 = d.inverse_cdf(p0).unwrap();
        let x1 = d.inverse_cdf(p1).unwrap();
        assert!(x0 <= 0.0 + 1e-6, "p=0 should give x near 0, got {x0}");
        assert!(x1 >= 1.0 - 1e-6, "p=1 should give x near 1, got {x1}");
    }

    #[test]
    fn inverse_cdf_boundary_probabilities_discrete() {
        let d = TenPoint;
        // p=0.0: cdf(0) = 0.1 >= 0.0, so leftmost satisfying k = 0
        let p0 = Probability::new(0.0).unwrap();
        assert_eq!(d.inverse_cdf(p0).unwrap(), 0u64);
        // p=1.0: cdf(9) = 1.0 >= 1.0, so leftmost satisfying k = 9
        let p1 = Probability::new(1.0).unwrap();
        assert_eq!(d.inverse_cdf(p1).unwrap(), 9u64);
    }
}
