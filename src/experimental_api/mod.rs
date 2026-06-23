//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature. All types here are subject to
//! breaking changes without a major version bump while in experimental status.

// mod bisect;  // Task 3
// mod traits;  // Task 4
mod types;

pub use types::{
    CdfError, InverseCdfError, InvalidDensity, InvalidMass, InvalidProbability, Probability,
    ProbabilityDensity, ProbabilityMass,
};
// pub use bisect::{bisection_search, Interval, PartitionSpace, SearchDirection};  // Task 3
// pub use traits::{Cdf, Pdf, Pmf};  // Task 4
