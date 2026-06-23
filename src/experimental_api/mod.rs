//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature. All types here are subject to
//! breaking changes without a major version bump while in experimental status.

pub(crate) mod bisect;
pub mod traits;
pub mod types;

pub use bisect::{Interval, PartitionSpace, SearchDirection};
pub use traits::{Cdf, DEFAULT_MAX_ITER, Pdf, Pmf};
pub use types::{
    CdfError, InverseCdfError, InvalidDensity, InvalidMass, InvalidProbability, Probability,
    ProbabilityDensity, ProbabilityMass,
};
