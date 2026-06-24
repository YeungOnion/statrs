pub(crate) trait PartitionSpace: Sized {
    type Point;
    type Cut;

    fn compute_cut(&self) -> Self::Cut;
    fn split(self, cut: &Self::Cut) -> (Self, Self);
    fn is_atomic(&self) -> bool;
    fn converge(self) -> Self::Point;
    fn point_from_cut(cut: Self::Cut) -> Self::Point;
}

pub(crate) enum SearchDirection {
    Left,
    Right,
    Found,
}

pub(crate) fn bisection_search<S>(
    mut space: S,
    oracle: impl Fn(&S::Cut) -> SearchDirection,
    max_iter: usize,
) -> Option<S::Point>
where
    S: PartitionSpace,
{
    for _ in 0..max_iter {
        if space.is_atomic() {
            return Some(space.converge());
        }
        let cut = space.compute_cut();
        match oracle(&cut) {
            SearchDirection::Found => return Some(S::point_from_cut(cut)),
            SearchDirection::Left => space = space.split(&cut).0,
            SearchDirection::Right => space = space.split(&cut).1,
        }
    }
    Some(space.converge())
}

#[derive(Copy, Clone, Debug)]
pub(crate) struct Interval<T> {
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
        (Interval { lo: self.lo, hi: *mid }, Interval { lo: *mid, hi: self.hi })
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
        (Interval { lo: self.lo, hi: *mid }, Interval { lo: *mid, hi: self.hi })
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

pub(crate) const DEFAULT_MAX_ITER: usize = 64;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisect_u64_finds_integer_sqrt() {
        let space = Interval { lo: 0u64, hi: 100u64 };
        let result = bisection_search(space, |cut| {
            match cut.checked_mul(*cut) {
                Some(sq) if sq == 25 => SearchDirection::Found,
                Some(sq) if sq > 25 => SearchDirection::Left,
                _ => SearchDirection::Right,
            }
        }, DEFAULT_MAX_ITER);
        assert_eq!(result, Some(5));
    }

    #[test]
    fn bisect_u64_atomic_terminates() {
        let space = Interval { lo: 3u64, hi: 4u64 };
        let result = bisection_search(space, |cut| {
            match cut.checked_mul(*cut) {
                Some(sq) if sq == 9 => SearchDirection::Found,
                Some(sq) if sq > 9 => SearchDirection::Left,
                _ => SearchDirection::Right,
            }
        }, DEFAULT_MAX_ITER);
        assert_eq!(result, Some(3));
    }

    #[test]
    fn bisect_f64_finds_threshold() {
        let space = Interval { lo: 0.0f64, hi: 1.0f64 };
        let result = bisection_search(space, |cut| {
            let diff = cut - 0.3;
            if diff.abs() < 1e-10 { SearchDirection::Found }
            else if diff > 0.0 { SearchDirection::Left }
            else { SearchDirection::Right }
        }, DEFAULT_MAX_ITER);
        assert!((result.unwrap() - 0.3).abs() < 1e-9);
    }

    #[test]
    fn interval_f64_is_never_atomic() {
        assert!(!Interval { lo: 0.0f64, hi: 1e-300f64 }.is_atomic());
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
