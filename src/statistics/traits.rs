use num_traits::float::Float;

/// The `Min` trait specifies than an object has a minimum value
pub trait Min<T> {
    /// Returns the minimum value in the domain of a given distribution
    /// if it exists, otherwise `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::statistics::Min;
    /// use statrs::distribution::{Uniform, UniformError};
    ///
    /// let n = Uniform::new(0.0, 1.0)?;
    /// assert_eq!(0.0, n.min());
    /// # Ok::<(), UniformError>(())
    /// ```
    fn min(&self) -> T;
}

/// The `Max` trait specifies that an object has a maximum value
pub trait Max<T> {
    /// Returns the maximum value in the domain of a given distribution
    /// if it exists, otherwise `None`.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::statistics::Max;
    /// use statrs::distribution::{Uniform, UniformError};
    ///
    /// let n = Uniform::new(0.0, 1.0)?;
    /// assert_eq!(1.0, n.max());
    /// # Ok::<(), UniformError>(())
    /// ```
    fn max(&self) -> T;
}

/// Exposes an entropy method, uses base e, (not Shannon, which base 2).
pub trait Entropy<T> {
    fn entropy(&self) -> T;
}

/// Trait to express co/variance as if it were a matrix,
///
/// For scalars this is variance and scalar
///
/// ```text
/// Sigma_ij = Cov[X_i, X_j]
/// ```
pub trait Covariance<T> {
    /// The dense form of this covariance matrix, doubly-indexed
    type M;
    /// The sparse/diagonal form of this covariance matrix, singly-indexed
    /// Can be described as potentially lossy form of the dense matrix
    type V;

    /// returns a covariance matrix, M_ij = Sigma_ij
    fn dense(&self) -> Self::M;

    /// returns a vector of marginal variances, v_i = Sigma_ii
    fn sparse(&self) -> Self::V;

    /// returns a vector scaled by covariance, (Sigma)^1/2 * vec{v}
    fn forward(&self, other: Self::V) -> Self::V;

    /// returns a vector unscaled by covariance, (Sigma)^-1/2 * vec{v}
    fn inverse(&self, other: Self::V) -> Self::V;

    /// returns the determinant of the covariance matrix, det(Sigma)
    fn determinant(&self) -> T;
}

#[cfg(feature = "nalgebra")]
mod multivariate {
    use nalgebra::{Cholesky, Dim, OMatrix, OVector, U1};

    impl<D> super::Covariance<f64> for OVector<f64, D>
    where
        D: Dim,
        nalgebra::DefaultAllocator:
            nalgebra::allocator::Allocator<f64, D> + nalgebra::allocator::Allocator<f64, D, D>,
    {
        type M = OMatrix<f64, D, D>;
        type V = Self;

        fn dense(&self) -> Self::M {
            OMatrix::from_diagonal(self)
        }

        fn sparse(&self) -> Self::V {
            self.clone()
        }

        fn forward(&self, other: Self::V) -> Self::V {
            let mut std_dev = self.clone().map(|x| x.sqrt());
            for i in 0..std_dev.len() {
                std_dev[i] *= other[i]
            }
            std_dev
        }

        fn inverse(&self, other: Self::V) -> Self::V {
            let mut inv_std_dev = self.clone().map(|x| x.sqrt().recip());
            for i in 0..inv_std_dev.len() {
                inv_std_dev[i] *= other[i]
            }
            inv_std_dev
        }

        fn determinant(&self) -> f64 {
            self.product()
        }
    }

    impl<D> super::Covariance<f64> for Cholesky<f64, D>
    where
        D: Dim,
        nalgebra::DefaultAllocator:
            nalgebra::allocator::Allocator<f64, D> + nalgebra::allocator::Allocator<f64, D, D>,
    {
        type M = OMatrix<f64, D, D>;
        type V = OVector<f64, D>;

        fn dense(&self) -> Self::M {
            self.l() * self.l().transpose()
        }

        fn sparse(&self) -> Self::V {
            let l = self.l();
            let (nr, _) = l.shape_generic();
            let mut res = OVector::zeros_generic(nr, U1);
            for (i, x) in res.iter_mut().enumerate() {
                *x = l[(i, i)];
            }
            res
        }

        fn forward(&self, other: Self::V) -> Self::V {
            self.l() * other
        }

        fn inverse(&self, other: Self::V) -> Self::V {
            self.solve(&other)
        }

        fn determinant(&self) -> f64 {
            self.determinant()
        }
    }
}

impl<T: Float> Covariance<T> for T {
    type M = T;
    type V = T;
    fn dense(&self) -> Self::M {
        *self
    }
    fn sparse(&self) -> Self::V {
        *self
    }
    fn forward(&self, other: Self::V) -> Self::V {
        self.sqrt() * other
    }
    fn inverse(&self, other: Self::V) -> Self::V {
        other / self.sqrt()
    }
    fn determinant(&self) -> T {
        *self
    }
}

/// Trait to express [standardized moments](https://en.wikipedia.org/wiki/Moment_(mathematics)#Standardized_moments)
/// of a distribution's realization implemented on its distribution's type
///
/// # Note for implementors
/// Associated types should capture semantics of the distribution, e.g.
/// [`Option`] should be used where moments is defined or undefined based on parameter values.
/// [`()`] is used for moments that are never defined
pub trait StandardizedMoment<T> {
    type Mu;
    type Var;
    type Kurt;
    type Skew;

    /// Returns the mean of a distribution, [defined as in measure theory](https://en.wikipedia.org/wiki/Expected_value#Arbitrary_real-valued_random_variables),
    ///
    /// # Formula
    /// ```text
    /// E[X], if defined
    /// ```
    /// In less specific terms, there are distributions with divergent integrals, which have no mean.
    fn mean(&self) -> Self::Mu;

    /// Returns the covariance of a distribution,
    ///
    ///
    /// # Formula
    /// ```text
    /// Cov[X_i, X_j] = E[(X - mu_i)(X_j - mu_j)]
    /// ```
    /// For univariate distributions this will be the variance.
    fn variance(&self) -> Self::Var;

    /// Returns the skewness,
    ///
    /// # Formula
    /// ```text
    /// E[((X - mu)/sigma)^3]
    /// ```
    fn skewness(&self) -> Self::Skew;

    /// Returns the excess kurtosis
    ///
    /// # Formula
    /// ```text
    /// E[((X - mu)/sigma)^4] - 3
    /// ```
    ///
    /// this is just the kurtosis minus 3
    fn excess_kurtosis(&self) -> Self::Kurt;
}

/// The `Mean` trait implements the calculation of a mean.
// TODO: Clarify the traits of multidimensional distributions
pub trait MeanN<T> {
    fn mean(&self) -> Option<T>;
}

// TODO: Clarify the traits of multidimensional distributions
pub trait VarianceN<T> {
    fn variance(&self) -> Option<T>;
}

/// The `Median` trait returns the median of the distribution.
pub trait Median<T> {
    /// Returns the median.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::statistics::Median;
    /// use statrs::distribution::{Uniform, UniformError};
    ///
    /// let n = Uniform::new(0.0, 1.0)?;
    /// assert_eq!(0.5, n.median());
    /// # Ok::<(), UniformError>(())
    /// ```
    fn median(&self) -> T;
}

/// The `Mode` trait specifies that an object has a closed form solution
/// for its mode(s)
pub trait Mode<T> {
    /// Returns the mode, if one exists.
    ///
    /// # Examples
    ///
    /// ```
    /// use statrs::statistics::Mode;
    /// use statrs::distribution::{Uniform, UniformError};
    ///
    /// let n = Uniform::new(0.0, 1.0)?;
    /// assert_eq!(Some(0.5), n.mode());
    /// # Ok::<(), UniformError>(())
    /// ```
    fn mode(&self) -> T;
}
