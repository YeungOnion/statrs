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
        (self.0.push(x), self.1.push(x), self.2.push(x), self.3.push(x))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Default)]
    struct Sum(f64);
    impl Accumulate for Sum {
        fn push(self, x: f64) -> Self { Sum(self.0 + x) }
    }

    #[derive(Default)]
    struct Count(u32);
    impl Accumulate for Count {
        fn push(self, _: f64) -> Self { Count(self.0 + 1) }
    }

    #[test]
    fn tuple2_fans_out() {
        let (s, c): (Sum, Count) = [1.0_f64, 2.0, 3.0]
            .iter().copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!(s.0, 6.0);
        assert_eq!(c.0, 3);
    }

    #[test]
    fn tuple3_fans_out() {
        let (a, b, c): (Sum, Sum, Sum) = [1.0_f64, 2.0]
            .iter().copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!(a.0, b.0);
        assert_eq!(b.0, c.0);
    }

    #[test]
    fn tuple4_fans_out() {
        let (a, b, c, d): (Sum, Sum, Sum, Sum) = [4.0_f64]
            .iter().copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!((a.0, b.0, c.0, d.0), (4.0, 4.0, 4.0, 4.0));
    }

    #[test]
    fn tuple5_fans_out() {
        let (a, b, c, d, e): (Sum, Sum, Sum, Sum, Sum) = [2.0_f64]
            .iter().copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!((a.0, b.0, c.0, d.0, e.0), (2.0, 2.0, 2.0, 2.0, 2.0));
    }

    #[test]
    fn single_tuple_fans_out() {
        let (s,): (Sum,) = [1.0_f64, 2.0]
            .iter().copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!(s.0, 3.0);
    }
}
