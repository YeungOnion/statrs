use super::{MultivariateNormal, density_distribution_exponential};
use crate::distribution::Continuous;
use crate::statistics::{Max, MeanN, Min, Mode, VarianceN};
use nalgebra::{Const, Dim, OMatrix, OVector};

impl<D> Min<OVector<f64, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn min(&self) -> OVector<f64, D> {
        OMatrix::repeat_generic(self.mu.shape_generic().0, Const::<1>, f64::NEG_INFINITY)
    }
}

impl<D> Max<OVector<f64, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn max(&self) -> OVector<f64, D> {
        OMatrix::repeat_generic(self.mu.shape_generic().0, Const::<1>, f64::INFINITY)
    }
}

impl<D> MeanN<OVector<f64, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn mean(&self) -> Option<OVector<f64, D>> {
        Some(self.mu.clone())
    }
}

impl<D> VarianceN<OMatrix<f64, D, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn variance(&self) -> Option<OMatrix<f64, D, D>> {
        Some(self.cov.clone())
    }
}

impl<D> Mode<OVector<f64, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn mode(&self) -> OVector<f64, D> {
        self.mu.clone()
    }
}

impl<D> Continuous<&OVector<f64, D>, f64> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn pdf(&self, x: &OVector<f64, D>) -> f64 {
        self.pdf_const
            * density_distribution_exponential(&self.mu, &self.precision, x)
                .unwrap()
                .exp()
    }

    fn ln_pdf(&self, x: &OVector<f64, D>) -> f64 {
        self.pdf_const.ln()
            + density_distribution_exponential(&self.mu, &self.precision, x).unwrap()
    }
}

#[rustfmt::skip]
#[cfg(test)]
mod tests  {
    use core::fmt::Debug;
    use crate::prec;
    use nalgebra::{dmatrix, dvector, matrix, vector, DimMin, OMatrix, OVector};

    use crate::{
        distribution::{Continuous, MultivariateNormal},
        statistics::{Max, MeanN, Min, Mode, VarianceN},
    };

    fn try_create<D>(mean: OVector<f64, D>, covariance: OMatrix<f64, D, D>) -> MultivariateNormal<D>
    where
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        let mvn = MultivariateNormal::new_from_nalgebra(mean, covariance);
        assert!(mvn.is_ok());
        mvn.unwrap()
    }

    fn create_case<D>(mean: OVector<f64, D>, covariance: OMatrix<f64, D, D>)
    where
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        let mvn = try_create(mean.clone(), covariance.clone());
        assert_eq!(mean, mvn.mean().unwrap());
        assert_eq!(covariance, mvn.variance().unwrap());
    }

    fn bad_create_case<D>(mean: OVector<f64, D>, covariance: OMatrix<f64, D, D>)
    where
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        let mvn = MultivariateNormal::new_from_nalgebra(mean, covariance);
        assert!(mvn.is_err());
    }

    fn test_case<T, F, D>(
        mean: OVector<f64, D>, covariance: OMatrix<f64, D, D>, expected: T, eval: F,
    ) where
        T: Debug + PartialEq,
        F: FnOnce(MultivariateNormal<D>) -> T,
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        let mvn = try_create(mean, covariance);
        let x = eval(mvn);
        assert_eq!(expected, x);
    }

    fn test_almost<F, D>(
        mean: OVector<f64, D>, covariance: OMatrix<f64, D, D>, expected: f64, acc: f64, eval: F,
    ) where
        F: FnOnce(MultivariateNormal<D>) -> f64,
        D: DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        let mvn = try_create(mean, covariance);
        let x = eval(mvn);
        prec::assert_abs_diff_eq!(expected, x, epsilon = acc);
    }

    #[test]
    fn test_create() {
        create_case(vector![0., 0.], matrix![1., 0.; 0., 1.]);
        create_case(vector![10., 5.], matrix![2., 1.; 1., 2.]);
        create_case(
            vector![4., 5., 6.],
            matrix![2., 1., 0.; 1., 2., 1.; 0., 1., 2.],
        );
        create_case(dvector![0., f64::INFINITY], dmatrix![1., 0.; 0., 1.]);
        create_case(
            dvector![0., 0.],
            dmatrix![f64::INFINITY, 0.; 0., f64::INFINITY],
        );
    }

    #[test]
    fn test_bad_create() {
        // Covariance not symmetric
        bad_create_case(vector![0., 0.], matrix![1., 1.; 0., 1.]);
        // Covariance not positive-definite
        bad_create_case(vector![0., 0.], matrix![1., 2.; 2., 1.]);
        // NaN in mean
        bad_create_case(dvector![0., f64::NAN], dmatrix![1., 0.; 0., 1.]);
        // NaN in Covariance Matrix
        bad_create_case(dvector![0., 0.], dmatrix![1., 0.; 0., f64::NAN]);
    }

    #[test]
    fn test_variance() {
        let variance = |x: MultivariateNormal<_>| x.variance().unwrap();
        test_case(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            matrix![1., 0.; 0., 1.],
            variance,
        );
        test_case(
            vector![0., 0.],
            matrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            matrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            variance,
        );
    }

    #[test]
    fn test_mode() {
        let mode = |x: MultivariateNormal<_>| x.mode();
        test_case(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            vector![0., 0.],
            mode,
        );
        test_case(
            vector![f64::INFINITY, f64::INFINITY],
            matrix![1., 0.; 0., 1.],
            vector![f64::INFINITY, f64::INFINITY],
            mode,
        );
    }

    #[test]
    fn test_min_max() {
        let min = |x: MultivariateNormal<_>| x.min();
        let max = |x: MultivariateNormal<_>| x.max();
        test_case(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            dvector![f64::NEG_INFINITY, f64::NEG_INFINITY],
            min,
        );
        test_case(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            dvector![f64::INFINITY, f64::INFINITY],
            max,
        );
        test_case(
            dvector![10., 1.],
            dmatrix![1., 0.; 0., 1.],
            dvector![f64::NEG_INFINITY, f64::NEG_INFINITY],
            min,
        );
        test_case(
            dvector![-3., 5.],
            dmatrix![1., 0.; 0., 1.],
            dvector![f64::INFINITY, f64::INFINITY],
            max,
        );
    }

    #[test]
    fn test_pdf() {
        let pdf = |arg| move |x: MultivariateNormal<_>| x.pdf(&arg);
        test_case(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            0.05854983152431917,
            pdf(vector![1., 1.]),
        );
        test_almost(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            0.013064233284684921,
            1e-15,
            pdf(vector![1., 2.]),
        );
        test_almost(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            1.8618676045881531e-23,
            1e-35,
            pdf(vector![1., 10.]),
        );
        test_almost(
            vector![0., 0.],
            matrix![1., 0.; 0., 1.],
            5.920684802611216e-45,
            1e-58,
            pdf(vector![10., 10.]),
        );
        test_almost(
            vector![0., 0.],
            matrix![1., 0.9; 0.9, 1.],
            1.6576716577547003e-05,
            1e-18,
            pdf(vector![1., -1.]),
        );
        test_almost(
            vector![0., 0.],
            matrix![1., 0.99; 0.99, 1.],
            4.1970621773477824e-44,
            1e-54,
            pdf(vector![1., -1.]),
        );
        test_almost(
            vector![0.5, -0.2],
            matrix![2.0, 0.3; 0.3, 0.5],
            0.0013075203140666656,
            1e-15,
            pdf(vector![2., 2.]),
        );
        test_case(
            vector![0., 0.],
            matrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            0.0,
            pdf(vector![10., 10.]),
        );
        test_case(
            vector![0., 0.],
            matrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            0.0,
            pdf(vector![100., 100.]),
        );
    }

    #[test]
    fn test_ln_pdf() {
        let ln_pdf = |arg| move |x: MultivariateNormal<_>| x.ln_pdf(&arg);
        test_case(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            (0.05854983152431917f64).ln(),
            ln_pdf(dvector![1., 1.]),
        );
        test_almost(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            (0.013064233284684921f64).ln(),
            1e-15,
            ln_pdf(dvector![1., 2.]),
        );
        test_almost(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            (1.8618676045881531e-23f64).ln(),
            1e-15,
            ln_pdf(dvector![1., 10.]),
        );
        test_almost(
            dvector![0., 0.],
            dmatrix![1., 0.; 0., 1.],
            (5.920684802611216e-45f64).ln(),
            1e-15,
            ln_pdf(dvector![10., 10.]),
        );
        test_almost(
            dvector![0., 0.],
            dmatrix![1., 0.9; 0.9, 1.],
            (1.6576716577547003e-05f64).ln(),
            1e-14,
            ln_pdf(dvector![1., -1.]),
        );
        test_almost(
            dvector![0., 0.],
            dmatrix![1., 0.99; 0.99, 1.],
            (4.1970621773477824e-44f64).ln(),
            1e-12,
            ln_pdf(dvector![1., -1.]),
        );
        test_almost(
            dvector![0.5, -0.2],
            dmatrix![2.0, 0.3; 0.3, 0.5],
            (0.0013075203140666656f64).ln(),
            1e-15,
            ln_pdf(dvector![2., 2.]),
        );
        test_case(
            dvector![0., 0.],
            dmatrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            f64::NEG_INFINITY,
            ln_pdf(dvector![10., 10.]),
        );
        test_case(
            dvector![0., 0.],
            dmatrix![f64::INFINITY, 0.; 0., f64::INFINITY],
            f64::NEG_INFINITY,
            ln_pdf(dvector![100., 100.]),
        );
    }

    #[test]
    #[should_panic]
    fn test_pdf_mismatched_arg_size() {
        let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.,]).unwrap();
        mvn.pdf(&vec![1.].into()); // x.size != mu.size
    }
}
