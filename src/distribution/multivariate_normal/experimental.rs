use super::{MultivariateNormal, density_distribution_exponential};
use crate::experimental_api::{
    ClosedFormCdf, InvalidVariate, Mean, Pdf, Probability, ProbabilityDensity, TryVariate,
    Variance, Variate,
};
use nalgebra::{DimMin, OVector};

impl<D> crate::Sealed for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
}

impl<D> TryVariate for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    type Repr = OVector<f64, D>;

    fn try_variate(
        &self,
        x: OVector<f64, D>,
    ) -> Result<Variate<Self, OVector<f64, D>>, InvalidVariate<OVector<f64, D>>> {
        if x.shape_generic().0 == self.mu.shape_generic().0 && x.iter().all(|v| v.is_finite()) {
            Ok(Variate::new(x))
        } else {
            Err(InvalidVariate(x))
        }
    }
}

impl<D> Pdf for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    fn pdf(&self, x: Variate<Self, OVector<f64, D>>) -> ProbabilityDensity {
        let exponential =
            density_distribution_exponential(&self.mu, &self.precision, &x.into_inner())
                .expect("dimensions already matched by try_variate");
        ProbabilityDensity::new(self.pdf_const * exponential.exp())
            .expect("multivariate normal PDF is always non-negative")
    }
}

impl<D> Mean for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    type Output = OVector<f64, D>;
    fn mean(&self) -> OVector<f64, D> {
        self.mu.clone()
    }
}

impl<D> Variance for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    type Output = OVector<f64, D>;
    fn variance(&self) -> OVector<f64, D> {
        self.cov.diagonal()
    }
}

impl<D> ClosedFormCdf for MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    /// No closed-form or reference implementation of the multivariate normal
    /// CDF exists in this crate yet.
    fn cdf(&self, _x: Variate<Self, OVector<f64, D>>) -> Probability {
        unimplemented!("multivariate normal CDF has no reference implementation in statrs yet")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_api::StdDev;
    use nalgebra::{dmatrix, dvector, matrix, vector};

    fn try_create<D>(
        mean: OVector<f64, D>,
        covariance: nalgebra::OMatrix<f64, D, D>,
    ) -> MultivariateNormal<D>
    where
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        MultivariateNormal::new_from_nalgebra(mean, covariance).unwrap()
    }

    #[test]
    fn try_variate_accepts_matching_dimension() {
        let d = try_create(vector![0., 0.], matrix![1., 0.; 0., 1.]);
        assert!(d.try_variate(vector![1., 1.]).is_ok());
    }

    #[test]
    fn try_variate_rejects_mismatched_dimension() {
        let d = try_create(dvector![0., 0.], dmatrix![1., 0.; 0., 1.]);
        assert!(d.try_variate(dvector![1.]).is_err());
    }

    #[test]
    fn try_variate_rejects_non_finite_component() {
        let d = try_create(vector![0., 0.], matrix![1., 0.; 0., 1.]);
        assert!(d.try_variate(vector![f64::NAN, 0.]).is_err());
    }

    #[test]
    fn pdf_matches_reference_value() {
        let d = try_create(vector![0., 0.], matrix![1., 0.; 0., 1.]);
        let x = d.try_variate(vector![1., 1.]).unwrap();
        crate::prec::assert_abs_diff_eq!(d.pdf(x).into_inner(), 0.05854983152431917, epsilon = 1e-15);
    }

    #[test]
    fn mean_matches_constructor_argument() {
        let d = try_create(vector![1., 2.], matrix![1., 0.; 0., 1.]);
        assert_eq!(d.mean(), vector![1., 2.]);
    }

    #[test]
    fn variance_is_covariance_diagonal_not_full_matrix() {
        let d = try_create(vector![0., 0.], matrix![2., 0.5; 0.5, 3.]);
        assert_eq!(d.variance(), vector![2., 3.]);
    }

    #[test]
    fn std_dev_is_elementwise_sqrt_of_variance() {
        let d = try_create(vector![0., 0.], matrix![4., 0.; 0., 9.]);
        assert_eq!(d.std_dev(), vector![2., 3.]);
    }

    #[test]
    #[should_panic]
    fn cdf_is_unimplemented() {
        let d = try_create(vector![0., 0.], matrix![1., 0.; 0., 1.]);
        let x = d.try_variate(vector![1., 1.]).unwrap();
        let _ = d.cdf(x);
    }
}
