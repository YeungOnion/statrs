//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature.

pub mod bisect;
pub mod traits;
pub mod types;

pub use bisect::{Interval, PartitionSpace, SearchDirection};
pub use traits::{Cdf, Pdf, Pmf};
pub use types::{
    CdfError, InvalidDensity, InvalidMass, InvalidProbability, InverseCdfError, Probability,
    ProbabilityDensity, ProbabilityMass,
};
