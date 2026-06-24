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
    /// Closed-form `inverse_cdf` override could not locate the target.
    NoConvergence,
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
        return Self(1.0 - self.0);
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
        }
    }
}
impl core::fmt::Display for InverseCdfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            InverseCdfError::NoConvergence => write!(f, "inverse_cdf bisection did not converge"),
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
