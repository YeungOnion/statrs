use core::f64;
use core::f64::consts::{E, PI};
#[cfg(feature = "rand")]
use nalgebra::Const;
use nalgebra::{Cholesky, DMatrix, DVector, Dim, DimMin, Dyn, OMatrix, OVector};

/// Computes both the normalization and exponential argument in the normal
/// distribution, returning `None` on dimension mismatch.
pub(super) fn density_normalization_and_exponential<D>(
    mu: &OVector<f64, D>,
    cov: &OMatrix<f64, D, D>,
    precision: &OMatrix<f64, D, D>,
    x: &OVector<f64, D>,
) -> Option<(f64, f64)>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    Some((
        density_distribution_pdf_const(mu, cov)?,
        density_distribution_exponential(mu, precision, x)?,
    ))
}

/// Computes the argument of the exponential term in the normal distribution,
/// returning `None` on dimension mismatch.
#[inline]
fn density_distribution_exponential<D>(
    mu: &OVector<f64, D>,
    precision: &OMatrix<f64, D, D>,
    x: &OVector<f64, D>,
) -> Option<f64>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    if x.shape_generic().0 != precision.shape_generic().0
        || x.shape_generic().0 != mu.shape_generic().0
        || !precision.is_square()
    {
        return None;
    }

    let dv = x - mu;
    let exp_term: f64 = -0.5 * (precision * &dv).dot(&dv);
    Some(exp_term)
}

/// Computes the argument of the normalization term in the normal distribution,
/// returning `None` on dimension mismatch.
#[inline]
fn density_distribution_pdf_const<D>(mu: &OVector<f64, D>, cov: &OMatrix<f64, D, D>) -> Option<f64>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    if cov.shape_generic().0 != mu.shape_generic().0 || !cov.is_square() {
        return None;
    }
    let cov_det = cov.determinant();
    Some(
        ((2. * PI).powi(mu.nrows() as i32) * cov_det.abs())
            .recip()
            .sqrt(),
    )
}

/// Implements the [Multivariate Normal](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
/// distribution using the "nalgebra" crate for matrix operations
///
/// # Examples
///
/// ```
/// use statrs::distribution::MultivariateNormal;
/// use nalgebra::{matrix, vector};
///
/// let mvn = MultivariateNormal::new_from_nalgebra(vector![0., 0.], matrix![1., 0.; 0., 1.]).unwrap();
///
/// #[cfg(feature = "experimental_api")]
/// {
///     use statrs::experimental_api::{Mean, Pdf, TryVariate};
///
///     assert_eq!(mvn.mean(), vector![0., 0.]);
///     let x = mvn.try_variate(vector![1., 1.]).unwrap();
///     assert_eq!(mvn.pdf(x).into_inner(), 0.05854983152431917);
/// }
///
/// #[cfg(not(feature = "experimental_api"))]
/// {
///     use statrs::distribution::Continuous;
///     use statrs::statistics::{MeanN, VarianceN};
///
///     assert_eq!(mvn.mean().unwrap(), vector![0., 0.]);
///     assert_eq!(mvn.variance().unwrap(), matrix![1., 0.; 0., 1.]);
///     assert_eq!(mvn.pdf(&vector![1., 1.]), 0.05854983152431917);
/// }
/// ```
#[derive(Clone, PartialEq, Debug)]
pub struct MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    cov_chol_decomp: OMatrix<f64, D, D>,
    mu: OVector<f64, D>,
    cov: OMatrix<f64, D, D>,
    precision: OMatrix<f64, D, D>,
    pdf_const: f64,
}

/// Represents the errors that can occur when creating a [`MultivariateNormal`].
#[derive(Copy, Clone, PartialEq, Eq, Debug, Hash)]
#[non_exhaustive]
pub enum MultivariateNormalError {
    /// The covariance matrix is asymmetric or contains a NaN.
    CovInvalid,

    /// The mean vector contains a NaN.
    MeanInvalid,

    /// The amount of rows in the vector of means is not equal to the amount
    /// of rows in the covariance matrix.
    DimensionMismatch,

    /// After all other validation, computing the Cholesky decomposition failed.
    CholeskyFailed,
}

impl core::fmt::Display for MultivariateNormalError {
    #[cfg_attr(coverage_nightly, coverage(off))]
    fn fmt(&self, f: &mut core::fmt::Formatter) -> core::fmt::Result {
        match self {
            MultivariateNormalError::CovInvalid => {
                write!(f, "Covariance matrix is asymmetric or contains a NaN")
            }
            MultivariateNormalError::MeanInvalid => write!(f, "Mean vector contains a NaN"),
            MultivariateNormalError::DimensionMismatch => write!(
                f,
                "Mean vector and covariance matrix do not have the same number of rows"
            ),
            MultivariateNormalError::CholeskyFailed => {
                write!(f, "Computing the Cholesky decomposition failed")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MultivariateNormalError {}

type NormalizedConstructorArguments<D> = (OVector<f64, D>, OMatrix<f64, D, D>, Cholesky<f64, D>);

impl MultivariateNormal<Dyn> {
    /// Constructs a new multivariate normal distribution with a mean of `mean`
    /// and covariance matrix `cov`
    ///
    /// # Errors
    ///
    /// Returns an error if the given covariance matrix is not
    /// symmetric or positive-definite
    pub fn new(mean: Vec<f64>, cov: Vec<f64>) -> Result<Self, MultivariateNormalError> {
        let mean = DVector::from_vec(mean);
        let cov = DMatrix::from_vec(mean.len(), mean.len(), cov);
        MultivariateNormal::new_from_nalgebra(mean, cov)
    }
}

/// Check that a covariance is square, perfectly symmetric, and non-NaN
fn check_cov<D>(cov: &OMatrix<f64, D, D>) -> Result<(), MultivariateNormalError>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    if !cov.is_square()
        || cov.lower_triangle() != cov.upper_triangle().transpose()
        || cov.iter().any(|f| f.is_nan())
    {
        Err(MultivariateNormalError::CovInvalid)
    } else {
        Ok(())
    }
}

/// Check the mean, covariance, and cholesky decomposition for incompatibilities, and return all three.
///
/// Covariance and cholesky decomposition are computed from each other as necessary.
///
/// # Panics
/// If both the `cov` and `cholesky` arguments are `None`; at least one must be `Some(_)`.
fn normalize_constructor_arguments<D>(
    mean: OVector<f64, D>,
    covariance: Option<OMatrix<f64, D, D>>,
    cholesky: Option<Cholesky<f64, D>>,
) -> Result<NormalizedConstructorArguments<D>, MultivariateNormalError>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    // Check that mean is valid
    if mean.iter().any(|f| f.is_nan()) {
        return Err(MultivariateNormalError::MeanInvalid);
    }
    let n = mean.shape_generic().0;

    let cov: OMatrix<f64, D, D>;
    let chol: Cholesky<f64, D>;

    match (covariance, cholesky) {
        (None, None) => {
            panic!("Must pass either cov or cholesky to normalize_constructor_arguments")
        }
        (Some(c), None) => {
            // Check covariance and compute Cholesky
            check_cov(&c)?;
            if c.shape_generic().0 != n {
                return Err(MultivariateNormalError::DimensionMismatch);
            }
            chol = Cholesky::new(c.clone()).ok_or(MultivariateNormalError::CholeskyFailed)?;
            cov = c;
        }
        (None, Some(ch)) => {
            // Check cholesky and compute covariance
            let ch_inner = ch.unpack_dirty();
            if ch_inner.shape_generic().0 != n {
                return Err(MultivariateNormalError::DimensionMismatch);
            }
            chol = Cholesky::pack_dirty(ch_inner);

            let l = chol.l();
            cov = l.clone() * l.transpose();
        }
        (Some(c), Some(ch)) => {
            // Check both covariance and cholesky
            check_cov(&c)?;
            if c.shape_generic().0 != n {
                return Err(MultivariateNormalError::DimensionMismatch);
            }
            cov = c;

            let ch_inner = ch.unpack_dirty();
            if ch_inner.shape_generic().0 != n {
                return Err(MultivariateNormalError::DimensionMismatch);
            }
            chol = Cholesky::pack_dirty(ch_inner);
        }
    }

    Ok((mean, cov, chol))
}

impl<D> MultivariateNormal<D>
where
    D: DimMin<D, Output = D>,
    nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
        + nalgebra::allocator::Allocator<D, D>
        + nalgebra::allocator::Allocator<D>,
{
    /// Constructs a new multivariate normal distribution with a mean of `mean`
    /// and covariance matrix `cov` using `nalgebra` `OVector` and `OMatrix`
    /// instead of `Vec<f64>`
    ///
    /// # Errors
    ///
    /// Returns an error if the given covariance matrix is not
    /// symmetric or positive-definite
    pub fn new_from_nalgebra(
        mean: OVector<f64, D>,
        cov: OMatrix<f64, D, D>,
    ) -> Result<Self, MultivariateNormalError> {
        let (mean, cov, cholesky) = normalize_constructor_arguments(mean, Some(cov), None)?;
        Ok(Self::new_unchecked(mean, cov, cholesky))
    }

    /// Construct a new multivariate normal from a mean and an already-computed
    /// Cholesky decomposition of the covariance matrix, using `nalgebra` types.
    ///
    /// Unlike [`MultivariateNormal::new_from_nalgebra`], this does not require
    /// the covariance matrix to be perfectly symmetric, since [`Cholesky`] is
    /// created with only the lower diagonal.
    ///
    /// # Errors
    ///
    /// Returns an error if the mean has any `NaN` values or the
    /// mean and cholesky decomposition have a different number of rows.
    pub fn new_from_cholesky(
        mean: OVector<f64, D>,
        cholesky: Cholesky<f64, D>,
    ) -> Result<Self, MultivariateNormalError> {
        let (mean, cov, cholesky) = normalize_constructor_arguments(mean, None, Some(cholesky))?;
        Ok(Self::new_unchecked(mean, cov, cholesky))
    }

    /// Construct a multivariate normal without checking the compatibility of the arguments.
    /// It is on the caller to ensure that they have the same shape and meet their respective
    /// invariants.
    fn new_unchecked(
        mean: OVector<f64, D>,
        cov: OMatrix<f64, D, D>,
        cholesky: Cholesky<f64, D>,
    ) -> MultivariateNormal<D> {
        // Grab precision
        let precision = cholesky.inverse();

        MultivariateNormal {
            pdf_const: density_distribution_pdf_const(&mean, &cov).unwrap(),
            cov_chol_decomp: cholesky.unpack(),
            mu: mean,
            cov,
            precision,
        }
    }

    /// Returns the entropy of the multivariate normal distribution
    ///
    /// # Formula
    ///
    /// ```text
    /// (1 / 2) * ln(det(2 * π * e * Σ))
    /// ```
    ///
    /// where `Σ` is the covariance matrix and `det` is the determinant
    pub fn entropy(&self) -> Option<f64> {
        Some(0.5 * self.cov.clone().scale(2. * PI * E).determinant().ln())
    }

    /// Returns the Cholesky decomposition of the covariance matrix
    ///
    /// # Example
    ///
    /// ```
    /// use nalgebra::OMatrix;
    /// use statrs::distribution::MultivariateNormal;
    ///
    /// let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();
    /// assert_eq!(mvn.clone_cov_chol_decomp().shape(), (2, 2));
    /// ```
    pub fn clone_cov_chol_decomp(&self) -> OMatrix<f64, D, D> {
        self.cov_chol_decomp.clone()
    }

    /// Returns the mean of the multivariate normal distribution
    ///
    /// # Example
    ///
    /// ```
    /// use nalgebra::{OVector, U2};
    /// use statrs::distribution::MultivariateNormal;
    ///
    /// let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();
    /// assert_eq!(mvn.mu(), &OVector::<f64, U2>::from_vec(vec![0., 0.]));
    /// ```
    pub fn mu(&self) -> &OVector<f64, D> {
        &self.mu
    }

    /// Returns the covariance of the multivariate normal distribution
    ///
    /// # Example
    ///
    /// ```
    /// use nalgebra::{OMatrix, U2};
    /// use statrs::distribution::MultivariateNormal;
    ///
    /// let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();
    /// assert_eq!(mvn.cov(), &OMatrix::<f64, U2, U2>::from_vec(vec![1., 0., 0., 1.]));
    /// ```
    pub fn cov(&self) -> &OMatrix<f64, D, D> {
        &self.cov
    }

    /// Returns the precision matrix of the multivariate normal distribution
    ///
    /// # Example
    ///
    /// ```
    /// use nalgebra::OMatrix;
    /// use statrs::distribution::MultivariateNormal;
    ///
    /// let mvn = MultivariateNormal::new(vec![0., 0.], vec![1., 0., 0., 1.]).unwrap();
    /// assert_eq!(mvn.precision().shape(), (2, 2));
    /// ```
    pub fn precision(&self) -> &OMatrix<f64, D, D> {
        &self.precision
    }
}

#[cfg(feature = "experimental_api")]
mod experimental;
#[cfg(not(feature = "experimental_api"))]
mod legacy;

impl<D> core::fmt::Display for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "N({}, {})", &self.mu, &self.cov)
    }
}

#[cfg(feature = "rand")]
#[cfg_attr(docsrs, doc(cfg(feature = "rand")))]
impl<D> ::rand::distr::Distribution<OVector<f64, D>> for MultivariateNormal<D>
where
    D: Dim,
    nalgebra::DefaultAllocator:
        nalgebra::allocator::Allocator<D> + nalgebra::allocator::Allocator<D, D>,
{
    /// Samples from the multivariate normal distribution
    ///
    /// # Formula
    /// ```text
    /// L * Z + μ
    /// ```
    ///
    /// where `L` is the Cholesky decomposition of the covariance matrix,
    /// `Z` is a vector of normally distributed random variables, and
    /// `μ` is the mean vector
    fn sample<R: ::rand::Rng + ?Sized>(&self, rng: &mut R) -> OVector<f64, D> {
        let d = crate::distribution::Normal::new(0., 1.).unwrap();
        let z = OVector::from_distribution_generic(self.mu.shape_generic().0, Const::<1>, &d, rng);
        (&self.cov_chol_decomp * z) + &self.mu
    }
}

#[cfg(test)]
mod tests {
    use crate::distribution::MultivariateNormal;
    use nalgebra::{dmatrix, dvector};

    fn try_create<D>(
        mean: nalgebra::OVector<f64, D>,
        covariance: nalgebra::OMatrix<f64, D, D>,
    ) -> MultivariateNormal<D>
    where
        D: nalgebra::DimMin<D, Output = D>,
        nalgebra::DefaultAllocator: nalgebra::allocator::Allocator<D>
            + nalgebra::allocator::Allocator<D, D>
            + nalgebra::allocator::Allocator<D>,
    {
        MultivariateNormal::new_from_nalgebra(mean, covariance).unwrap()
    }

    #[test]
    fn test_entropy() {
        let entropy = |x: MultivariateNormal<_>| x.entropy().unwrap();
        assert_eq!(
            entropy(try_create(dvector![0., 0.], dmatrix![1., 0.; 0., 1.])),
            2.8378770664093453
        );
        assert_eq!(
            entropy(try_create(dvector![0., 0.], dmatrix![1., 0.5; 0.5, 1.])),
            2.694036030183455
        );
        assert_eq!(
            entropy(try_create(
                dvector![0., 0.],
                dmatrix![f64::INFINITY, 0.; 0., f64::INFINITY]
            )),
            f64::INFINITY
        );
    }

    #[test]
    fn test_error_is_sync_send() {
        fn assert_sync_send<T: Sync + Send>() {}
        assert_sync_send::<super::MultivariateNormalError>();
    }
}
