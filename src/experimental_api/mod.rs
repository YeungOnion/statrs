//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature.
//! Expect this to change rapidly, potentially breaking every commit.

pub mod bisect;
pub mod traits;
pub mod types;

pub use traits::{ClosedFormCdf, InverseCdf, Pdf, Pmf, TryVariate};
pub use types::{
    InvalidDensity, InvalidMass, InvalidProbability, InvalidVariate, InverseCdfError, Probability,
    ProbabilityDensity, ProbabilityMass, Variate,
};
