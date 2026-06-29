/// A type that accumulates `f64` observations one at a time.
///
/// Tuple impls up to arity 5 fan each observation out to all accumulators,
/// so multiple statistics can share a single iterator pass:
///
/// ```
/// use statrs::experimental_api::{Accumulate, Moments, VarianceAccum, AbsMinAccum};
///
/// let data = [3.0_f64, -1.0, 4.0, 1.0, -5.0];
/// let (var, mn): (VarianceAccum, AbsMinAccum) = data.iter().copied()
///     .fold(Default::default(), Accumulate::push);
///
/// assert_eq!(mn.abs_min(), Some(1.0));
/// assert!((var.variance().unwrap() - 12.8).abs() < 1e-10);
/// ```
pub trait Accumulate: Default + Sized {
    fn push(self, x: f64) -> Self;
}

impl<A: Accumulate> Accumulate for (A,) {
    fn push(self, x: f64) -> Self {
        (self.0.push(x),)
    }
}

impl<A: Accumulate, B: Accumulate> Accumulate for (A, B) {
    fn push(self, x: f64) -> Self {
        (self.0.push(x), self.1.push(x))
    }
}

impl<A: Accumulate, B: Accumulate, C: Accumulate> Accumulate for (A, B, C) {
    fn push(self, x: f64) -> Self {
        (self.0.push(x), self.1.push(x), self.2.push(x))
    }
}

impl<A: Accumulate, B: Accumulate, C: Accumulate, D: Accumulate> Accumulate for (A, B, C, D) {
    fn push(self, x: f64) -> Self {
        (
            self.0.push(x),
            self.1.push(x),
            self.2.push(x),
            self.3.push(x),
        )
    }
}

impl<A: Accumulate, B: Accumulate, C: Accumulate, D: Accumulate, E: Accumulate> Accumulate
    for (A, B, C, D, E)
{
    fn push(self, x: f64) -> Self {
        (
            self.0.push(x),
            self.1.push(x),
            self.2.push(x),
            self.3.push(x),
            self.4.push(x),
        )
    }
}
