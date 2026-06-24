//! Validated newtypes for probability values and distribution domains.
//!
//! [`Domain<D>`] ties a validated `f64` to a specific distribution via a phantom type
//! parameter. Validate once at the boundary with `TryFrom`; after that the compiler
//! enforces both that raw `f64` can't slip through and that values from one distribution
//! can't be passed to another's methods.
//!
//! ```compile_fail
//! # use statrs::experimental_api::{Cdf, Domain, HasSupport, Probability};
//! # struct A;
//! # impl HasSupport for A { type Bound = f64; fn contains(x: f64) -> bool { x.is_finite() } }
//! # impl Cdf for A { fn cdf(&self, x: Domain<Self>) -> Probability { todo!() } }
//! A.cdf(0.5_f64); // expected Domain<A>, found f64
//! ```
//!
//! ```compile_fail
//! # use statrs::experimental_api::{Cdf, Domain, HasSupport, Probability};
//! # struct A;
//! # impl HasSupport for A { type Bound = f64; fn contains(x: f64) -> bool { x.is_finite() } }
//! # impl Cdf for A { fn cdf(&self, x: Domain<Self>) -> Probability { todo!() } }
//! # struct B;
//! # impl HasSupport for B { type Bound = f64; fn contains(x: f64) -> bool { x.is_finite() } }
//! let x: Domain<B> = 0.5_f64.try_into().unwrap();
//! A.cdf(x); // Domain<B> is not Domain<A>
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

/// Error returned by CDF computations.
#[non_exhaustive]
pub enum CdfError {
    InvalidInput,
}

/// Error returned by inverse-CDF computations.
#[non_exhaustive]
pub enum InverseCdfError {
    /// Used for search failures where inverse_cdf is not closed form.
    NoConvergence,
    /// The exact inverse maps to a point outside the distribution's support
    /// (e.g. `p = 1.0` for a half-open upper bound).
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

impl core::fmt::Debug for CdfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CdfError::InvalidInput => write!(f, "CdfError::InvalidInput"),
        }
    }
}
impl core::fmt::Display for CdfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CdfError::InvalidInput => write!(f, "invalid input to cdf"),
        }
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
impl std::error::Error for CdfError {}
#[cfg(feature = "std")]
impl std::error::Error for InverseCdfError {}

// ---- Domain validation ----

/// Static support membership for distributions whose bounds are type-level constants.
///
/// `Bound` is the type of values that can assert membership and act as bisection
/// boundaries. All near-term impls use `Bound = f64`.
pub trait HasSupport {
    type Bound: Copy;
    fn contains(x: Self::Bound) -> bool;
}

/// A value known to lie within distribution `D`'s support.
///
/// Construct via `TryFrom<D::Bound>` or the equivalent `try_into()`.
pub struct Domain<D: HasSupport>(D::Bound, PhantomData<D>);

impl<D: HasSupport> Domain<D> {
    pub fn into_inner(self) -> D::Bound {
        self.0
    }
}

impl<D: HasSupport> Copy for Domain<D> {}
impl<D: HasSupport> Clone for Domain<D> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<D: HasSupport> PartialEq for Domain<D>
where
    D::Bound: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}
impl<D: HasSupport> core::fmt::Debug for Domain<D>
where
    D::Bound: core::fmt::Debug,
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Domain").field(&self.0).finish()
    }
}

/// Error returned when a value falls outside a distribution's support.
pub struct InvalidDomain<B>(pub B);

impl<B: core::fmt::Display> core::fmt::Debug for InvalidDomain<B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "InvalidDomain({})", self.0)
    }
}
impl<B: core::fmt::Display> core::fmt::Display for InvalidDomain<B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "value {} is outside the distribution's support", self.0)
    }
}
#[cfg(feature = "std")]
impl<B: core::fmt::Display> std::error::Error for InvalidDomain<B> {}

impl<D: HasSupport<Bound = f64>> TryFrom<f64> for Domain<D> {
    type Error = InvalidDomain<f64>;
    fn try_from(x: f64) -> Result<Self, Self::Error> {
        if D::contains(x) {
            Ok(Self(x, PhantomData))
        } else {
            Err(InvalidDomain(x))
        }
    }
}

impl<D: HasSupport<Bound = u64>> TryFrom<u64> for Domain<D> {
    type Error = InvalidDomain<u64>;
    fn try_from(x: u64) -> Result<Self, Self::Error> {
        if D::contains(x) {
            Ok(Self(x, PhantomData))
        } else {
            Err(InvalidDomain(x))
        }
    }
}

impl TryFrom<f64> for Probability {
    type Error = InvalidProbability;
    fn try_from(v: f64) -> Result<Self, Self::Error> {
        Self::new(v)
    }
}


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
