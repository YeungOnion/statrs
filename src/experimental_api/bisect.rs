pub trait PartitionSpace: Sized {
    type Point;
    type Cut;

    fn compute_cut(&self) -> Self::Cut;
    fn split(self, cut: &Self::Cut) -> (Self, Self);

    /// Returns true when the region cannot be meaningfully subdivided.
    /// Checked before `compute_cut` each iteration to avoid infinite loops
    /// on discrete domains where `lo + 1 == hi`.
    fn is_atomic(&self) -> bool;

    /// Extracts a representative point from a narrowed region.
    fn converge(self) -> Self::Point;

    /// Converts a cut value into the `Point` type.
    fn point_from_cut(cut: Self::Cut) -> Self::Point;
}

pub enum SearchDirection {
    Left,
    Right,
    Found,
}

pub(crate) trait SearchOracle<C> {
    fn evaluate(&self, cut: &C) -> SearchDirection;
}

pub(crate) fn bisection_search<S, O>(
    mut space: S,
    oracle: &O,
    max_iterations: usize,
) -> Option<S::Point>
where
    S: PartitionSpace,
    O: SearchOracle<S::Cut>,
{
    for _ in 0..max_iterations {
        if space.is_atomic() {
            return Some(space.converge());
        }
        let cut = space.compute_cut();
        match oracle.evaluate(&cut) {
            SearchDirection::Found => return Some(S::point_from_cut(cut)),
            SearchDirection::Left => space = space.split(&cut).0,
            SearchDirection::Right => space = space.split(&cut).1,
        }
    }
    Some(space.converge())
}

#[derive(Copy, Clone, Debug)]
pub struct Interval<T> {
    pub lo: T,
    pub hi: T,
}

impl PartitionSpace for Interval<f64> {
    type Point = f64;
    type Cut = f64;

    fn compute_cut(&self) -> f64 {
        self.lo + (self.hi - self.lo) / 2.0
    }

    fn split(self, mid: &f64) -> (Self, Self) {
        (
            Interval {
                lo: self.lo,
                hi: *mid,
            },
            Interval {
                lo: *mid,
                hi: self.hi,
            },
        )
    }

    fn is_atomic(&self) -> bool {
        false
    }

    fn converge(self) -> f64 {
        self.lo + (self.hi - self.lo) / 2.0
    }

    fn point_from_cut(cut: f64) -> f64 {
        cut
    }
}

impl PartitionSpace for Interval<u64> {
    type Point = u64;
    type Cut = u64;

    fn compute_cut(&self) -> u64 {
        self.lo + (self.hi - self.lo) / 2
    }

    fn split(self, mid: &u64) -> (Self, Self) {
        (
            Interval {
                lo: self.lo,
                hi: *mid,
            },
            Interval {
                lo: *mid,
                hi: self.hi,
            },
        )
    }

    fn is_atomic(&self) -> bool {
        self.hi.saturating_sub(self.lo) <= 1
    }

    fn converge(self) -> u64 {
        self.lo
    }

    fn point_from_cut(cut: u64) -> u64 {
        cut
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Oracle that finds the integer square root: smallest k with k*k >= n
    struct IntSqrtOracle {
        n: u64,
    }
    impl SearchOracle<u64> for IntSqrtOracle {
        fn evaluate(&self, cut: &u64) -> SearchDirection {
            match cut.checked_mul(*cut) {
                Some(sq) if sq == self.n => SearchDirection::Found,
                Some(sq) if sq > self.n => SearchDirection::Left,
                _ => SearchDirection::Right,
            }
        }
    }

    #[test]
    fn bisect_u64_finds_integer_sqrt() {
        let space = Interval { lo: 0u64, hi: 100u64 };
        let oracle = IntSqrtOracle { n: 25 };
        let result = bisection_search(space, &oracle, 64);
        assert_eq!(result, Some(5));
    }

    #[test]
    fn bisect_u64_atomic_terminates() {
        let space = Interval { lo: 3u64, hi: 4u64 };
        let oracle = IntSqrtOracle { n: 9 };
        let result = bisection_search(space, &oracle, 64);
        assert_eq!(result, Some(3));
    }

    // Oracle that finds where a monotone f64 function crosses a threshold
    struct ThresholdOracle {
        threshold: f64,
    }
    impl SearchOracle<f64> for ThresholdOracle {
        fn evaluate(&self, cut: &f64) -> SearchDirection {
            let diff = cut - self.threshold;
            if diff.abs() < 1e-10 {
                SearchDirection::Found
            } else if diff > 0.0 {
                SearchDirection::Left
            } else {
                SearchDirection::Right
            }
        }
    }

    #[test]
    fn bisect_f64_finds_threshold() {
        let space = Interval { lo: 0.0f64, hi: 1.0f64 };
        let oracle = ThresholdOracle { threshold: 0.3 };
        let result = bisection_search(space, &oracle, 64).unwrap();
        assert!((result - 0.3).abs() < 1e-9);
    }

    #[test]
    fn interval_f64_is_never_atomic() {
        let space = Interval { lo: 0.0f64, hi: 1e-300f64 };
        assert!(!space.is_atomic());
    }

    #[test]
    fn interval_u64_atomic_when_width_one() {
        assert!(Interval { lo: 5u64, hi: 6u64 }.is_atomic());
        assert!(!Interval { lo: 5u64, hi: 7u64 }.is_atomic());
    }

    #[test]
    fn interval_u64_atomic_when_equal() {
        assert!(Interval { lo: 5u64, hi: 5u64 }.is_atomic());
    }
}

use core::marker::PhantomData;
use crate::experimental_api::types::Probability;
use crate::experimental_api::traits::Cdf;

pub(crate) struct CdfOracle<'a, D, K> {
    dist: &'a D,
    target: Probability,
    epsilon: f64,
    _k: PhantomData<K>,
}

impl<'a, D, K> CdfOracle<'a, D, K> {
    pub(crate) fn new(dist: &'a D, target: Probability) -> Self {
        Self { dist, target, epsilon: 1e-10, _k: PhantomData }
    }
}

impl<D: Cdf<f64>> SearchOracle<f64> for CdfOracle<'_, D, f64> {
    fn evaluate(&self, cut: &f64) -> SearchDirection {
        match self.dist.cdf(*cut) {
            Err(_) => SearchDirection::Right,
            Ok(p) => {
                let diff = p.into_inner() - self.target.into_inner();
                if diff.abs() < self.epsilon {
                    SearchDirection::Found
                } else if diff > 0.0 {
                    SearchDirection::Left
                } else {
                    SearchDirection::Right
                }
            }
        }
    }
}

impl<D: Cdf<u64>> SearchOracle<u64> for CdfOracle<'_, D, u64> {
    fn evaluate(&self, cut: &u64) -> SearchDirection {
        match self.dist.cdf(*cut) {
            Err(_) => SearchDirection::Right,
            Ok(p) if p >= self.target => {
                // Check if this is the leftmost satisfying point.
                // If cut == 0 or cdf(cut - 1) < target, we have found the minimum.
                let prev_satisfies = cut
                    .checked_sub(1)
                    .and_then(|prev| self.dist.cdf(prev).ok())
                    .map(|pp| pp >= self.target)
                    .unwrap_or(false);
                if prev_satisfies {
                    SearchDirection::Left
                } else {
                    SearchDirection::Found
                }
            }
            Ok(_) => SearchDirection::Right,
        }
    }
}
