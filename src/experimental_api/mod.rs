//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature.
//! Expect this to change rapidly, potentially breaking every commit.

pub mod bisect;
pub mod traits;
pub mod types;

pub use bisect::{Interval, PartitionSpace, SearchDirection};
pub use traits::{Cdf, InverseCdf, Pdf, Pmf};
pub use types::{
    CdfError, Domain, HasSupport, InvalidDensity, InvalidDomain, InvalidMass, InvalidProbability,
    InverseCdfError, Probability, ProbabilityDensity, ProbabilityMass,
};
