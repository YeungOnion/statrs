use crate::Sealed;

mod private_sqrt {
    pub(crate) trait Sqrt {
        fn sqrt(self) -> Self;
    }
    impl Sqrt for f64 {
        fn sqrt(self) -> Self {
            self.sqrt()
        }
    }
    impl<F, const N: usize> Sqrt for [F; N]
    where
        F: Sqrt,
    {
        fn sqrt(self) -> Self {
            self.map(Sqrt::sqrt)
        }
    }
    impl<T: Sqrt> Sqrt for Option<T> {
        fn sqrt(self) -> Self {
            self.map(Sqrt::sqrt)
        }
    }
}
use private_sqrt::Sqrt;

/// The mean of a distribution or accumulator.
///
/// `Output` is picked by the implementor: `f64` when the mean always exists,
/// `Option<f64>` when it exists only under some condition on the parameters
/// or amount of data seen. A type whose mean never exists (e.g. Cauchy)
/// simply does not implement this trait.
#[allow(private_bounds)]
pub trait Mean: Sealed {
    type Output;
    fn mean(&self) -> Self::Output;
}

/// Sample variance, normalised by N-1 (Bessel's correction).
///
/// See [`PopulationVariance`] for the N-normalised form. `Output` follows the
/// same convention as [`Mean::Output`].
#[allow(private_bounds)]
pub trait Variance: Sealed {
    type Output;
    fn variance(&self) -> Self::Output;
}

/// Standard deviation, derived from [`Variance::variance`].
///
/// Implemented automatically for any type whose `Variance::Output` is built
/// out of `f64` (elementwise, however deeply nested in arrays/`Option`) — see
/// the blanket impl below. No manual impl is needed or allowed.
#[allow(private_bounds)]
pub trait StdDev: Sealed {
    type Output;
    fn std_dev(&self) -> Self::Output;
}

impl<T> StdDev for T
where
    T: Variance,
    T::Output: Sqrt,
{
    type Output = T::Output;
    fn std_dev(&self) -> Self::Output {
        self.variance().sqrt()
    }
}

/// Population variance, normalised by N rather than N-1.
///
/// Use when the data is the full population. For sample estimates see
/// [`Variance::variance`].
#[allow(private_bounds)]
pub trait PopulationVariance: Sealed {
    /// Population variance, normalised by N.
    fn population_variance(&self) -> Option<f64>;
    fn population_std_dev(&self) -> Option<f64> {
        self.population_variance().map(f64::sqrt)
    }
}

#[allow(private_bounds)]
pub trait Skewness: Sealed {
    fn skewness(&self) -> Option<f64>;
}

#[allow(private_bounds)]
pub trait Entropy: Sealed {
    fn entropy(&self) -> Option<f64>;
}

#[allow(private_bounds)]
pub trait Median: Sealed {
    fn median(&self) -> Option<f64>;
}

#[allow(private_bounds)]
pub trait Mode: Sealed {
    fn mode(&self) -> Option<f64>;
}

#[allow(private_bounds)]
pub trait Min: Sealed {
    fn min(&self) -> Option<f64>;
}

#[allow(private_bounds)]
pub trait Max: Sealed {
    fn max(&self) -> Option<f64>;
}

#[cfg(test)]
mod tests {
    use super::*;

    struct ConstDist {
        mean: f64,
        variance: f64,
    }

    impl Sealed for ConstDist {}

    impl Mean for ConstDist {
        type Output = Option<f64>;
        fn mean(&self) -> Self::Output {
            Some(self.mean)
        }
    }

    impl Variance for ConstDist {
        type Output = Option<f64>;
        fn variance(&self) -> Self::Output {
            Some(self.variance)
        }
    }

    struct NoVariance;

    impl Sealed for NoVariance {}

    impl Mean for NoVariance {
        type Output = Option<f64>;
        fn mean(&self) -> Self::Output {
            Some(1.0)
        }
    }

    impl Variance for NoVariance {
        type Output = Option<f64>;
        fn variance(&self) -> Self::Output {
            None
        }
    }

    struct ScalarDist;

    impl Sealed for ScalarDist {}

    impl Variance for ScalarDist {
        type Output = f64;
        fn variance(&self) -> Self::Output {
            9.0
        }
    }

    #[test]
    fn std_dev_derived_from_variance() {
        let d = ConstDist {
            mean: 0.0,
            variance: 4.0,
        };
        assert_eq!(d.std_dev(), Some(2.0));
    }

    #[test]
    fn std_dev_none_when_variance_none() {
        assert_eq!(NoVariance.std_dev(), None);
    }

    #[test]
    fn std_dev_scalar_output() {
        assert_eq!(ScalarDist.std_dev(), 3.0);
    }

    #[test]
    fn min_impl_returns_some() {
        struct HasMin;
        impl Sealed for HasMin {}
        impl Min for HasMin {
            fn min(&self) -> Option<f64> {
                Some(1.0)
            }
        }
        assert_eq!(HasMin.min(), Some(1.0));
    }

    #[test]
    fn max_impl_returns_none_for_empty() {
        struct EmptyMax;
        impl Sealed for EmptyMax {}
        impl Max for EmptyMax {
            fn max(&self) -> Option<f64> {
                None
            }
        }
        assert_eq!(EmptyMax.max(), None);
    }
}
