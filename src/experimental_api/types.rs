//! Validated newtypes for probability values and distribution variates.
//!
//! Call [`TryVariate::try_variate`] on a distribution to get a [`Variate<D, B>`].
//! The phantom `D` prevents it from reaching a different distribution's
//! methods; raw values can't reach [`ClosedFormCdf::cdf`], [`Pdf::pdf`], etc. directly.
//!
//! [`Variate`] construction is crate-internal; external [`TryVariate`]
//! implementations are not yet supported.
//!
//! ```compile_fail
//! # use statrs::experimental_api::{ClosedFormCdf, Variate, Probability, TryVariate, InvalidVariate};
//! # struct A;
//! # impl TryVariate for A {
//! #     type Bound = f64;
//! #     fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> { todo!() }
//! # }
//! # impl ClosedFormCdf for A { fn cdf(&self, x: Variate<Self, f64>) -> Probability { todo!() } }
//! let a = A;
//! a.cdf(0.5_f64); // expected Variate<A, f64>, found f64
//! ```
//!
//! ```compile_fail
//! # use statrs::experimental_api::{ClosedFormCdf, Variate, Probability, TryVariate, InvalidVariate};
//! # struct A;
//! # impl TryVariate for A {
//! #     type Bound = f64;
//! #     fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> { todo!() }
//! # }
//! # impl ClosedFormCdf for A { fn cdf(&self, x: Variate<Self, f64>) -> Probability { todo!() } }
//! # struct B;
//! # impl TryVariate for B {
//! #     type Bound = f64;
//! #     fn try_variate(&self, x: f64) -> Result<Variate<Self, f64>, InvalidVariate<f64>> { todo!() }
//! # }
//! let b = B;
//! let x: Variate<B, f64> = b.try_variate(0.5).unwrap();
//! let a = A;
//! a.cdf(x); // Variate<B, f64> is not Variate<A, f64>
//! ```

use core::marker::PhantomData;
use decorum::R64;

/// A probability value in `[0.0, 1.0]`, guaranteed finite.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Probability(R64);

/// A probability density value in `[0.0, ∞)`, guaranteed finite.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ProbabilityDensity(R64);

/// A probability mass value in `[0.0, 1.0]`, guaranteed finite.
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ProbabilityMass(R64);

/// Error returned when a value is not a valid probability.
pub struct InvalidProbability(pub f64);

/// Error returned when a value is not a valid probability density.
pub struct InvalidDensity(pub f64);

/// Error returned when a value is not a valid probability mass.
pub struct InvalidMass(pub f64);

/// Error returned by inverse-CDF computations.
#[non_exhaustive]
pub enum InverseCdfError {
    /// The bisection search did not converge within the iteration limit.
    NoConvergence,
    /// The inverse falls outside the support, e.g. `p = 1.0` with a
    /// half-open upper bound.
    OutOfSupport,
}

impl Probability {
    pub fn new(v: f64) -> Result<Self, InvalidProbability> {
        let f = R64::try_new(v).map_err(|_| InvalidProbability(v))?;
        if !(0.0..=1.0).contains(&f.into_inner()) {
            return Err(InvalidProbability(v));
        }
        Ok(Self(f))
    }

    pub fn complement(&self) -> Self {
        Self(1.0 - self.0)
    }

    pub fn into_inner(self) -> f64 {
        self.0.into_inner()
    }
}

impl ProbabilityDensity {
    pub fn new(v: f64) -> Result<Self, InvalidDensity> {
        let f = R64::try_new(v).map_err(|_| InvalidDensity(v))?;
        if f.into_inner() < 0.0 {
            return Err(InvalidDensity(v));
        }
        Ok(Self(f))
    }

    pub fn into_inner(self) -> f64 {
        self.0.into_inner()
    }
}

impl ProbabilityMass {
    pub fn new(v: f64) -> Result<Self, InvalidMass> {
        let f = R64::try_new(v).map_err(|_| InvalidMass(v))?;
        if !(0.0..=1.0).contains(&f.into_inner()) {
            return Err(InvalidMass(v));
        }
        Ok(Self(f))
    }

    pub fn into_inner(self) -> f64 {
        self.0.into_inner()
    }
}

impl core::fmt::Debug for InvalidProbability {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InvalidProbability({})", self.0)
    }
}
impl core::fmt::Display for InvalidProbability {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "invalid probability: {}", self.0)
    }
}

impl core::fmt::Debug for InvalidDensity {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InvalidDensity({})", self.0)
    }
}
impl core::fmt::Display for InvalidDensity {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "invalid probability density: {}", self.0)
    }
}

impl core::fmt::Debug for InvalidMass {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InvalidMass({})", self.0)
    }
}
impl core::fmt::Display for InvalidMass {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "invalid probability mass: {}", self.0)
    }
}

impl core::fmt::Debug for InverseCdfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            InverseCdfError::NoConvergence => write!(f, "InverseCdfError::NoConvergence"),
            InverseCdfError::OutOfSupport => write!(f, "InverseCdfError::OutOfSupport"),
        }
    }
}
impl core::fmt::Display for InverseCdfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            InverseCdfError::NoConvergence => write!(f, "inverse_cdf search did not converge"),
            InverseCdfError::OutOfSupport => {
                write!(f, "inverse_cdf result falls outside the support")
            }
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for InvalidProbability {}
#[cfg(feature = "std")]
impl std::error::Error for InvalidDensity {}
#[cfg(feature = "std")]
impl std::error::Error for InvalidMass {}
#[cfg(feature = "std")]
impl std::error::Error for InverseCdfError {}

impl TryFrom<f64> for Probability {
    type Error = InvalidProbability;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        Self::new(v)
    }
}

// ---- Variate ----

/// A value of a random variable, validated against distribution `D`'s support.
///
/// `D` is a phantom type tying this value to a specific distribution —
/// `Variate<Beta, f64>` and `Variate<Normal, f64>` are distinct types even
/// though both wrap an `f64`.
///
/// `B` is the type backing the sample space:
/// - `f64` for continuous univariate
/// - `u64` for discrete
///
/// Construct via [`TryVariate::try_variate`].
#[doc(alias = "domain")]
#[doc(alias = "sample")]
pub struct Variate<D, B>(B, PhantomData<D>);

impl<D, B: Copy> Variate<D, B> {
    /// Returns the inner value, consuming the wrapper.
    pub fn into_inner(self) -> B {
        self.0
    }

    /// Constructs without a membership check; callers must ensure `x` is in-support.
    pub(crate) fn new(x: B) -> Self {
        Self(x, PhantomData)
    }
}

impl<D, B: Copy> Copy for Variate<D, B> {}
impl<D, B: Copy> Clone for Variate<D, B> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D, B: PartialEq> PartialEq for Variate<D, B> {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<D, B: core::fmt::Debug> core::fmt::Debug for Variate<D, B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Variate").field(&self.0).finish()
    }
}

/// Error returned when a value falls outside a distribution's support.
pub struct InvalidVariate<B>(pub B);

impl<B: core::fmt::Display> core::fmt::Debug for InvalidVariate<B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InvalidVariate({})", self.0)
    }
}
impl<B: core::fmt::Display> core::fmt::Display for InvalidVariate<B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "value {} is outside the distribution's support", self.0)
    }
}
#[cfg(feature = "std")]
impl<B: core::fmt::Display> std::error::Error for InvalidVariate<B> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probability_accepts_zero_and_one() {
        assert!(Probability::new(0.0).is_ok());
        assert!(Probability::new(1.0).is_ok());
    }

    #[test]
    fn probability_accepts_interior() {
        assert!(Probability::new(0.5).is_ok());
    }

    #[test]
    fn probability_rejects_negative() {
        assert!(Probability::new(-0.001).is_err());
    }

    #[test]
    fn probability_rejects_above_one() {
        assert!(Probability::new(1.001).is_err());
    }

    #[test]
    fn probability_rejects_nan() {
        assert!(Probability::new(f64::NAN).is_err());
    }

    #[test]
    fn probability_rejects_inf() {
        assert!(Probability::new(f64::INFINITY).is_err());
        assert!(Probability::new(f64::NEG_INFINITY).is_err());
    }

    #[test]
    fn probability_into_inner_roundtrips() {
        let p = Probability::new(0.75).unwrap();
        assert_eq!(p.into_inner(), 0.75);
    }

    #[test]
    fn density_accepts_zero_and_large() {
        assert!(ProbabilityDensity::new(0.0).is_ok());
        assert!(ProbabilityDensity::new(100.0).is_ok()); // density can exceed 1
    }

    #[test]
    fn density_rejects_negative() {
        assert!(ProbabilityDensity::new(-0.001).is_err());
    }

    #[test]
    fn density_rejects_nan_and_inf() {
        assert!(ProbabilityDensity::new(f64::NAN).is_err());
        assert!(ProbabilityDensity::new(f64::INFINITY).is_err());
    }

    #[test]
    fn mass_accepts_unit_interval() {
        assert!(ProbabilityMass::new(0.0).is_ok());
        assert!(ProbabilityMass::new(1.0).is_ok());
        assert!(ProbabilityMass::new(0.3).is_ok());
    }

    #[test]
    fn mass_rejects_out_of_range() {
        assert!(ProbabilityMass::new(-0.001).is_err());
        assert!(ProbabilityMass::new(1.001).is_err());
    }
}
