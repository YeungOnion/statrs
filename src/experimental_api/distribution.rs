use crate::Sealed;

/// Sample mean and variance. `std_dev` defaults to `sqrt(variance)`.
///
/// `variance` is normalised by N-1 (Bessel's correction). See
/// [`PopulationMoments`] for the N-normalised form.
#[allow(private_bounds)]
pub trait Moments: Sealed {
    fn mean(&self) -> Option<f64>;
    /// Sample variance, normalised by N-1.
    fn variance(&self) -> Option<f64>;
    fn std_dev(&self) -> Option<f64> {
        self.variance().map(f64::sqrt)
    }
}

/// Population variance, normalised by N rather than N-1.
///
/// See also [`Moments::variance`] for sample estimates.
#[allow(private_bounds)]
pub trait PopulationMoments: Sealed {
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

    impl Moments for ConstDist {
        fn mean(&self) -> Option<f64> {
            Some(self.mean)
        }
        fn variance(&self) -> Option<f64> {
            Some(self.variance)
        }
    }

    struct NoVariance;

    impl Sealed for NoVariance {}

    impl Moments for NoVariance {
        fn mean(&self) -> Option<f64> {
            Some(1.0)
        }
        fn variance(&self) -> Option<f64> {
            None
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
