//! Traits for CDF, PDF, PMF, inverse CDF, and variate validation.
//!
//! Methods take [`Variate<Self, Self::Repr>`] instead of raw values.
//! Validate once at the entry point with [`TryVariate::try_variate`]; the
//! methods are infallible from there.
//!
//! ```ignore
//! // Implementing TryVariate requires Variate::new, which is pub(crate).
//! // See the tests module in this file for a full worked example.
//! # use statrs::experimental_api::{ClosedFormCdf, InverseCdf, Variate, TryVariate, Probability, InvalidVariate};
//! # struct MyDist;
//! # impl TryVariate for MyDist {
//! #     fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
//! #         if x.is_finite() && (0.0..=1.0).contains(&x) {
//! #             Ok(Variate::new(x))
//! #         } else {
//! #             Err(InvalidVariate(x))
//! #         }
//! #     }
//! # }
//! # impl ClosedFormCdf for MyDist {
//! #     fn cdf(&self, x: Variate<Self, f64>) -> Probability {
//! #         Probability::new(x.into_inner()).unwrap()
//! #     }
//! # }
//! # impl InverseCdf for MyDist {
//! #     fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, f64>, InverseCdfError> {
//! #         self.try_variate(p.into_inner()).map_err(|_| InverseCdfError::OutOfSupport)
//! #     }
//! # }
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let d = MyDist;
//! let x: Variate<_> = d.try_variate(0.6)?;  // validate at the boundary
//! let p = d.cdf(x);                          // infallible from here
//! let x2 = d.inverse_cdf(p)?;               // Variate<MyDist, f64> — passes directly back to cdf
//! let _ = d.cdf(x2);
//! # Ok(())
//! # }
//! ```

use crate::experimental_api::types::{InvalidVariate, InverseCdfError, Probability, Variate};
use crate::experimental_api::{ProbabilityDensity, ProbabilityMass};

/// Validates a raw value against a distribution's support, returning a [`Variate`].
///
/// `Repr` is the backing primitive of the sample space:
/// - `f64` for continuous univariate continuous
/// - `u64` for discrete.
///
/// Call [`try_variate`][TryVariate::try_variate] at the boundary; pass the
/// returned [`Variate`] to [`ClosedFormCdf::cdf`], [`Pdf::pdf`], etc.
/// Passing it to a different distribution's methods is a compile error.
pub trait TryVariate: Sized {
    type Repr: Clone;
    fn try_variate(
        &self,
        x: Self::Repr,
    ) -> Result<Variate<Self, Self::Repr>, InvalidVariate<Self::Repr>>;
}

/// Probability density function.
///
/// `ln_pdf` defaults to `pdf(x).into_inner().ln()`.
pub trait Pdf: TryVariate {
    fn pdf(&self, x: Variate<Self, Self::Repr>) -> ProbabilityDensity;

    fn ln_pdf(&self, x: Variate<Self, Self::Repr>) -> f64 {
        self.pdf(x).into_inner().ln()
    }
}

/// Probability mass function.
///
/// `ln_pmf` defaults to `pmf(x).into_inner().ln()`.
pub trait Pmf: TryVariate {
    fn pmf(&self, x: Variate<Self, Self::Repr>) -> ProbabilityMass;

    fn ln_pmf(&self, x: Variate<Self, Self::Repr>) -> f64 {
        self.pmf(x).into_inner().ln()
    }
}

/// Cumulative distribution function.
///
/// Takes [`Variate<Self, Self::Repr>`] rather than a raw value, so the
/// method is infallible.
///
/// Note to implementors:
/// If the method for Cdf evaluation is fallible, that should be a new trait.
pub trait ClosedFormCdf: TryVariate {
    fn cdf(&self, x: Variate<Self, Self::Repr>) -> Probability;
}

/// Inverse CDF. Implement directly — closed form, or a distribution-specific search.
pub trait InverseCdf: ClosedFormCdf {
    fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, Self::Repr>, InverseCdfError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_api::types::{InvalidVariate, ProbabilityDensity, ProbabilityMass};

    // Continuous uniform on [0, 1].
    struct UnitUniform;

    impl TryVariate for UnitUniform {
        type Repr = f64;
        fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
            if x.is_finite() && (0.0..=1.0).contains(&x) {
                Ok(Variate::new(x))
            } else {
                Err(InvalidVariate(x))
            }
        }
    }

    impl Pdf for UnitUniform {
        fn pdf(&self, _x: Variate<Self, f64>) -> ProbabilityDensity {
            ProbabilityDensity::new(1.0).unwrap()
        }
    }

    // Discrete uniform on {0, 1, ..., 9}.
    struct TenPoint;

    impl TryVariate for TenPoint {
        type Repr = f64;
        fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
            if x.is_finite() && x >= 0.0 && x <= 9.0 && x.fract() == 0.0 {
                Ok(Variate::new(x))
            } else {
                Err(InvalidVariate(x))
            }
        }
    }

    impl Pmf for TenPoint {
        fn pmf(&self, _x: Variate<Self, f64>) -> ProbabilityMass {
            ProbabilityMass::new(0.1).unwrap()
        }
    }

    // Uniform on [1, 2) — used to test closed-form inverse_cdf.
    struct Uniform12;

    impl TryVariate for Uniform12 {
        type Repr = f64;
        fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> {
            if x.is_finite() && x >= 1.0 && x < 2.0 {
                Ok(Variate::new(x))
            } else {
                Err(InvalidVariate(x))
            }
        }
    }

    impl ClosedFormCdf for Uniform12 {
        fn cdf(&self, x: Variate<Self, f64>) -> Probability {
            Probability::new(x.into_inner() - 1.0).expect("x - 1 ∈ [0,1) for x ∈ [1,2)")
        }
    }

    impl InverseCdf for Uniform12 {
        fn inverse_cdf(&self, p: Probability) -> Result<Variate<Self, f64>, InverseCdfError> {
            self.try_variate(p.into_inner() + 1.0)
                .map_err(|_| InverseCdfError::OutOfSupport)
        }
    }

    #[test]
    fn pdf_unit_uniform_is_one() {
        let d = UnitUniform;
        let x = d.try_variate(0.5).unwrap();
        assert_eq!(d.pdf(x).into_inner(), 1.0);
    }

    #[test]
    fn ln_pdf_default_impl() {
        let d = UnitUniform;
        let x = d.try_variate(0.5).unwrap();
        assert!((d.ln_pdf(x) - 0.0_f64).abs() < 1e-12);
    }

    #[test]
    fn pdf_rejects_out_of_support() {
        assert!(UnitUniform.try_variate(1.5).is_err());
    }

    #[test]
    fn pmf_ten_point_is_one_tenth() {
        let d = TenPoint;
        let x = d.try_variate(3.0).unwrap();
        assert!((d.pmf(x).into_inner() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn ln_pmf_default_impl() {
        let d = TenPoint;
        let x = d.try_variate(7.0).unwrap();
        assert!((d.ln_pmf(x) - 0.1_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn pmf_rejects_non_integer() {
        assert!(TenPoint.try_variate(2.5).is_err());
    }

    #[test]
    fn cdf_midpoint() {
        let d = Uniform12;
        let x = d.try_variate(1.5).unwrap();
        assert_eq!(d.cdf(x).into_inner(), 0.5);
    }

    #[test]
    fn cdf_lower_boundary() {
        let x = Uniform12.try_variate(1.0).unwrap();
        assert_eq!(Uniform12.cdf(x).into_inner(), 0.0);
    }

    #[test]
    fn inverse_cdf_closed_form_roundtrip() {
        let d = Uniform12;
        let p = Probability::new(0.3).unwrap();
        let x = d.inverse_cdf(p).unwrap();
        assert!((d.cdf(x).into_inner() - 0.3).abs() < 1e-12);
    }

    #[test]
    fn inverse_cdf_rejects_p1() {
        let p = Probability::new(1.0).unwrap();
        assert!(Uniform12.inverse_cdf(p).is_err());
    }

    #[test]
    fn domain_rejects_out_of_support() {
        assert!(Uniform12.try_variate(0.5).is_err());
        assert!(Uniform12.try_variate(2.0).is_err());
        assert!(Uniform12.try_variate(f64::NAN).is_err());
    }

    #[test]
    fn domain_accepts_lower_rejects_upper_bound() {
        assert!(Uniform12.try_variate(1.0).is_ok());
        assert!(Uniform12.try_variate(2.0).is_err());
    }
}
