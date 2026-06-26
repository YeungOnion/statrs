use core::ops::ControlFlow;

pub trait FoldStat: Iterator<Item = f64> + Sized {
    fn fold_stat<S, O, E>(
        mut self,
        init: S,
        mut step: impl FnMut(S, f64) -> ControlFlow<Result<O, E>, S>,
        finalize: impl FnOnce(S) -> Option<O>,
    ) -> Result<Option<O>, E> {
        match self.try_fold(init, |s, x| step(s, x)) {
            ControlFlow::Continue(s) => Ok(finalize(s)),
            ControlFlow::Break(r) => r.map(Some),
        }
    }
}

impl<I: Iterator<Item = f64>> FoldStat for I {}

#[cfg(test)]
mod tests {
    use super::*;
    use core::ops::ControlFlow;

    fn add_step(acc: f64, x: f64) -> ControlFlow<Result<f64, f64>, f64> {
        if x < 0.0 {
            ControlFlow::Break(Err(x))
        } else {
            ControlFlow::Continue(acc + x)
        }
    }

    fn sum_finalize(acc: f64) -> Option<f64> {
        if acc == 0.0 { None } else { Some(acc) }
    }

    #[test]
    fn empty_iterator_returns_ok_none() {
        let result = []
            .iter()
            .copied()
            .fold_stat(0.0_f64, add_step, sum_finalize);
        assert_eq!(result, Ok(None));
    }

    #[test]
    fn valid_inputs_finalize() {
        let result = [1.0_f64, 2.0, 3.0]
            .iter()
            .copied()
            .fold_stat(0.0_f64, add_step, sum_finalize);
        assert_eq!(result, Ok(Some(6.0)));
    }

    #[test]
    fn domain_error_short_circuits() {
        let result =
            [1.0_f64, -2.0, 3.0]
                .iter()
                .copied()
                .fold_stat(0.0_f64, add_step, sum_finalize);
        assert_eq!(result, Err(-2.0));
    }
}
