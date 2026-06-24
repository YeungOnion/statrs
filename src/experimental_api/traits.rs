//! Trait definitions for CDF, PDF, PMF, and inverse CDF.
//!
//! All methods take [`Domain<Self>`][crate::experimental_api::Domain] rather than raw `f64`.
//! Validation happens once at the call boundary via `TryFrom`; from there `?` propagates
//! errors and the trait methods themselves are infallible or return only convergence errors.
//!
//! ```
//! # use statrs::experimental_api::{Cdf, InverseCdf, Domain, HasSupport, Probability};
//! # struct MyDist;
//! # impl HasSupport for MyDist {
//! #     type Bound = f64;
//! #     fn contains(x: f64) -> bool { x.is_finite() && (0.0..=1.0).contains(&x) }
//! # }
//! # impl Cdf for MyDist {
//! #     fn cdf(&self, x: Domain<Self>) -> Probability { Probability::new(x.into_inner()).unwrap() }
//! # }
//! # impl InverseCdf for MyDist {
//! #     fn search_bounds(&self) -> (f64, f64) { (0.0, 1.0) }
//! # }
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let d = MyDist;
//! let x: Domain<_> = 0.6_f64.try_into()?; // validate at the boundary
//! let p = d.cdf(x);                         // infallible
//! let x2 = d.inverse_cdf(p)?;              // Domain<MyDist> — passes directly back to cdf
//! let _ = d.cdf(x2);
//! # Ok(())
//! # }
//! ```
//!
//! The existing distributions delegate to their current implementations, so both APIs
//! agree on values. Here `Beta` shows the contrast — the old CDF takes raw `f64` and
//! can panic; the new one takes a validated `Domain<Beta>` and is infallible:
//!
//! ```
//! # use statrs::distribution::{Beta, ContinuousCDF};
//! # use statrs::experimental_api::{Cdf, Domain};
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let dist = Beta::new(2.0, 5.0)?;
//! let raw = 0.3_f64;
//!
//! let old: f64 = ContinuousCDF::cdf(&dist, raw);       // accepts f64, can panic
//! let x: Domain<Beta> = raw.try_into()?;
//! let new: f64 = Cdf::cdf(&dist, x).into_inner();      // accepts Domain<Beta>, infallible
//!
//! assert!((old - new).abs() < f64::EPSILON);
//! # Ok(())
//! # }
//! ```

use crate::experimental_api::bisect::{
    bisection_search, Interval, SearchDirection, DEFAULT_MAX_ITER,
};
use crate::experimental_api::types::{
    Domain, HasSupport, InverseCdfError, Probability, ProbabilityDensity, ProbabilityMass,
};

/// PDF for distributions whose support is encoded in the type.
///
/// `ln_pdf` has a default implementation via `pdf`.
pub trait Pdf: HasSupport + Sized {
    fn pdf(&self, x: Domain<Self>) -> ProbabilityDensity;

    fn ln_pdf(&self, x: Domain<Self>) -> f64 {
        self.pdf(x).into_inner().ln()
    }
}

/// PMF for distributions whose support is encoded in the type.
///
/// `ln_pmf` has a default implementation via `pmf`.
pub trait Pmf: HasSupport + Sized {
    fn pmf(&self, x: Domain<Self>) -> ProbabilityMass;

    fn ln_pmf(&self, x: Domain<Self>) -> f64 {
        self.pmf(x).into_inner().ln()
    }
}

/// CDF for distributions whose support is encoded in the type.
///
/// The return is infallible: `Domain<Self>` is proof the input is in-support.
pub trait Cdf: HasSupport + Sized {
    fn cdf(&self, x: Domain<Self>) -> Probability;
}

/// Inverse CDF extending [`Cdf`].
///
/// Implement `search_bounds` to get a bisection-based `inverse_cdf` by default,
/// or override it with a closed-form solution. Distributions whose search space
/// is not a scalar interval (e.g. Dirichlet) should override `inverse_cdf` directly
/// and leave `search_bounds` unreachable.
pub trait InverseCdf: Cdf {
    /// The `(lo, hi)` bounds within which bisection searches for the inverse.
    ///
    /// Must be a subset of the support. For half-open upper bounds (e.g. `[1, 2)`)
    /// use the open endpoint as `hi` — bisection never lands exactly there, and
    /// `Domain::try_from` rejects it if it does.
    fn search_bounds(&self) -> (f64, f64);

    fn inverse_cdf(&self, p: Probability) -> Result<Domain<Self>, InverseCdfError>
    where
        Self: HasSupport<Bound = f64>,
    {
        let (lo, hi) = self.search_bounds();
        bisection_search(Interval { lo, hi }, |cut| {
            let Ok(x) = Domain::try_from(*cut) else {
                return SearchDirection::Right;
            };
            let diff = self.cdf(x).into_inner() - p.into_inner();
            if diff.abs() < 1e-10 {
                SearchDirection::Found
            } else if diff > 0.0 {
                SearchDirection::Left
            } else {
                SearchDirection::Right
            }
        }, DEFAULT_MAX_ITER)
        .and_then(|x| Domain::try_from(x).ok())
        .ok_or(InverseCdfError::NoConvergence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::experimental_api::types::{
        InvalidDomain, Probability, ProbabilityDensity, ProbabilityMass,
    };

    struct Uniform12;

    impl HasSupport for Uniform12 {
        type Bound = f64;
        fn contains(x: f64) -> bool {
            x.is_finite() && x >= 1.0 && x < 2.0
        }
    }

    impl Cdf for Uniform12 {
        fn cdf(&self, x: Domain<Self>) -> Probability {
            Probability::new(x.into_inner() - 1.0).expect("x - 1 ∈ [0,1) for x ∈ [1,2)")
        }
    }

    impl InverseCdf for Uniform12 {
        fn search_bounds(&self) -> (f64, f64) {
            (1.0, 2.0)
        }

        fn inverse_cdf(&self, p: Probability) -> Result<Domain<Self>, InverseCdfError> {
            (p.into_inner() + 1.0)
                .try_into()
                .map_err(|_: InvalidDomain<f64>| InverseCdfError::OutOfSupport)
        }
    }

    // Continuous uniform on [0, 1].
    struct UnitUniform;
    impl HasSupport for UnitUniform {
        type Bound = f64;
        fn contains(x: f64) -> bool {
            x.is_finite() && (0.0..=1.0).contains(&x)
        }
    }
    impl Pdf for UnitUniform {
        fn pdf(&self, _x: Domain<Self>) -> ProbabilityDensity {
            ProbabilityDensity::new(1.0).unwrap()
        }
    }

    // Discrete uniform on {0, 1, ..., 9}.
    struct TenPoint;
    impl HasSupport for TenPoint {
        type Bound = f64;
        fn contains(x: f64) -> bool {
            x.is_finite() && x >= 0.0 && x <= 9.0 && x.fract() == 0.0
        }
    }
    impl Pmf for TenPoint {
        fn pmf(&self, _x: Domain<Self>) -> ProbabilityMass {
            ProbabilityMass::new(0.1).unwrap()
        }
    }

    #[test]
    fn pdf_unit_uniform_is_one() {
        let d = UnitUniform;
        let x: Domain<_> = 0.5_f64.try_into().unwrap();
        assert_eq!(d.pdf(x).into_inner(), 1.0);
    }

    #[test]
    fn ln_pdf_default_impl() {
        let d = UnitUniform;
        let x: Domain<_> = 0.5_f64.try_into().unwrap();
        assert!((d.ln_pdf(x) - 0.0_f64).abs() < 1e-12);
    }

    #[test]
    fn pdf_rejects_out_of_support() {
        let r: Result<Domain<UnitUniform>, _> = 1.5_f64.try_into();
        assert!(r.is_err());
    }

    #[test]
    fn pmf_ten_point_is_one_tenth() {
        let d = TenPoint;
        let x: Domain<_> = 3.0_f64.try_into().unwrap();
        assert!((d.pmf(x).into_inner() - 0.1).abs() < 1e-12);
    }

    #[test]
    fn ln_pmf_default_impl() {
        let d = TenPoint;
        let x: Domain<_> = 7.0_f64.try_into().unwrap();
        assert!((d.ln_pmf(x) - 0.1_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn pmf_rejects_non_integer() {
        let r: Result<Domain<TenPoint>, _> = 2.5_f64.try_into();
        assert!(r.is_err());
    }

    #[test]
    fn cdf_midpoint() {
        let d = Uniform12;
        let x: Domain<_> = 1.5_f64.try_into().unwrap();
        assert_eq!(d.cdf(x).into_inner(), 0.5);
    }

    #[test]
    fn cdf_lower_boundary() {
        let x: Domain<Uniform12> = 1.0_f64.try_into().unwrap();
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
    fn inverse_cdf_bisection_roundtrip() {
        struct Uniform12Bisect;
        impl HasSupport for Uniform12Bisect {
            type Bound = f64;
            fn contains(x: f64) -> bool {
                x.is_finite() && x >= 1.0 && x < 2.0
            }
        }
        impl Cdf for Uniform12Bisect {
            fn cdf(&self, x: Domain<Self>) -> Probability {
                Probability::new(x.into_inner() - 1.0).unwrap()
            }
        }
        impl InverseCdf for Uniform12Bisect {
            fn search_bounds(&self) -> (f64, f64) {
                (1.0, 2.0)
            }
        }

        let d = Uniform12Bisect;
        let p = Probability::new(0.3).unwrap();
        let x = d.inverse_cdf(p).unwrap();
        assert!((d.cdf(x).into_inner() - 0.3).abs() < 1e-6);
    }

    #[test]
    fn inverse_cdf_rejects_p1() {
        let p = Probability::new(1.0).unwrap();
        assert!(Uniform12.inverse_cdf(p).is_err());
    }

    #[test]
    fn domain_rejects_out_of_support() {
        let r: Result<Domain<Uniform12>, _> = 0.5_f64.try_into();
        assert!(r.is_err());
        let r: Result<Domain<Uniform12>, _> = 2.0_f64.try_into();
        assert!(r.is_err());
        let r: Result<Domain<Uniform12>, _> = f64::NAN.try_into();
        assert!(r.is_err());
    }

    #[test]
    fn domain_accepts_lower_rejects_upper_bound() {
        let lo: Result<Domain<Uniform12>, _> = 1.0_f64.try_into();
        assert!(lo.is_ok());
        let hi: Result<Domain<Uniform12>, _> = 2.0_f64.try_into();
        assert!(hi.is_err());
    }
}
