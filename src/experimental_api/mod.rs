//! Experimental non-panicking API for statistical distributions.
//!
//! Gated on the `experimental_api` feature.
//! Expect this to change rapidly, potentially breaking every commit.

pub mod bisect;
pub mod distribution;
pub mod fold;
pub mod streaming;
pub mod summary;
pub mod traits;
pub mod types;

pub use distribution::{Entropy, Max, Median, Min, Mode, Moments, PopulationMoments, Skewness};
pub use fold::Accumulate;
pub use streaming::{
    AbsMaxAccum, AbsMinAccum, CovAccum, MeanAccum, RunningCov, SkewnessAccum, VarianceAccum, kahan,
};
pub use summary::{abs_max, abs_min, covariance, geometric_mean, harmonic_mean, quadratic_mean};

pub use traits::{ClosedFormCdf, InverseCdf, Pdf, Pmf, TryVariate};
pub use types::{
    InvalidDensity, InvalidMass, InvalidProbability, InvalidVariate, InverseCdfError, Probability,
    ProbabilityDensity, ProbabilityMass, Variate,
};
