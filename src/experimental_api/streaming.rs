use crate::experimental_api::{Max, Min, Moments, Skewness};

pub struct RunningMoments<const ORDER: usize> {
    pub count: u64,
    min: f64,
    max: f64,
    m: [f64; ORDER],
}

impl<const ORDER: usize> Default for RunningMoments<ORDER> {
    fn default() -> Self {
        Self {
            count: 0,
            min: f64::INFINITY,
            max: f64::NEG_INFINITY,
            m: [0.0; ORDER],
        }
    }
}

impl<const ORDER: usize> RunningMoments<ORDER> {
    /// Folds one observation into the accumulator (Welford online algorithm).
    ///
    /// Designed for use as an `Iterator::fold` accumulator:
    /// ```
    /// # use statrs::experimental_api::VarianceAccum;
    /// let s = [1.0_f64, 2.0, 3.0].iter().copied()
    ///     .fold(VarianceAccum::default(), VarianceAccum::push);
    /// ```
    pub fn push(mut self, x: f64) -> Self {
        self.count += 1;
        let n = self.count as f64;

        // NaN in the data propagates outward; f64::min/max don't propagate NaN,
        // so we handle it explicitly.
        if x.is_nan() {
            self.min = f64::NAN;
            self.max = f64::NAN;
        } else if !self.min.is_nan() {
            if x < self.min { self.min = x; }
            if x > self.max { self.max = x; }
        }

        // Welford / Pébaÿ (2008) central moment update.
        // m[0] = running mean (M1)
        // m[1] = sum of squared deviations (M2), if ORDER >= 2
        // m[2] = sum of cubed deviations (M3), if ORDER >= 3
        //
        // Update order: M3 before M2 (uses old M2), M2 before M1.
        let delta = x - self.m[0];    // deviation from OLD mean
        let delta_n = delta / n;
        let new_mean = self.m[0] + delta_n;
        let delta2 = x - new_mean;    // deviation from NEW mean

        if let Some(old_m2) = self.m.get(1).copied() {
            if let Some(m3) = self.m.get_mut(2) {
                let delta_n2 = delta_n * delta_n;
                *m3 += delta * delta_n2 * (n - 1.0) * (n - 2.0) / n
                    - 3.0 * delta_n * old_m2;
            }
            if let Some(m2) = self.m.get_mut(1) {
                *m2 += delta * delta2;
            }
        }

        self.m[0] = new_mean;
        self
    }
}

impl crate::experimental_api::distribution::private::Sealed for RunningMoments<2> {}
impl crate::experimental_api::distribution::private::Sealed for RunningMoments<3> {}

impl Moments for RunningMoments<2> {
    fn mean(&self) -> Option<f64> {
        if self.count == 0 { None } else { Some(self.m[0]) }
    }

    fn variance(&self) -> Option<f64> {
        if self.count < 2 {
            None
        } else {
            Some(self.m[1] / (self.count - 1) as f64)
        }
    }
}

impl Min for RunningMoments<2> {
    fn min(&self) -> Option<f64> {
        if self.count == 0 { None } else { Some(self.min) }
    }
}

impl Max for RunningMoments<2> {
    fn max(&self) -> Option<f64> {
        if self.count == 0 { None } else { Some(self.max) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type Accum = RunningMoments<2>;

    #[test]
    fn default_is_empty() {
        let s = Accum::default();
        assert_eq!(s.count, 0);
        assert!(s.min.is_infinite() && s.min > 0.0);
        assert!(s.max.is_infinite() && s.max < 0.0);
    }

    #[test]
    fn push_increments_count() {
        let s = Accum::default().push(3.0);
        assert_eq!(s.count, 1);
    }

    #[test]
    fn push_single_element_min_max() {
        let s = Accum::default().push(7.0);
        assert_eq!(s.min, 7.0);
        assert_eq!(s.max, 7.0);
    }

    #[test]
    fn push_nan_propagates_to_min_max() {
        let s = Accum::default().push(1.0).push(f64::NAN).push(2.0);
        assert!(s.min.is_nan());
        assert!(s.max.is_nan());
    }

    #[test]
    fn push_two_elements_mean_accumulator() {
        // [2, 8] → mean = 5.0; m[0] should be 5.0
        let s = Accum::default().push(2.0).push(8.0);
        assert!((s.m[0] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn moments_empty_returns_none() {
        let s = RunningMoments::<2>::default();
        assert_eq!(s.mean(), None);
        assert_eq!(s.variance(), None);
        assert_eq!(s.std_dev(), None);
    }

    #[test]
    fn moments_single_element_variance_none() {
        let s = RunningMoments::<2>::default().push(5.0);
        assert_eq!(s.mean(), Some(5.0));
        assert_eq!(s.variance(), None);
        assert_eq!(s.std_dev(), None);
    }

    #[test]
    fn moments_known_dataset() {
        // [2,4,4,4,5,5,7,9]: mean=5.0, sample variance=32/7
        let data = [2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s = data.iter().copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!((s.mean().unwrap() - 5.0).abs() < 1e-12);
        assert!((s.variance().unwrap() - 32.0 / 7.0).abs() < 1e-12);
        assert!((s.std_dev().unwrap() - (32.0_f64 / 7.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn moments_nan_propagates_as_some_nan() {
        let s = [1.0_f64, f64::NAN]
            .iter().copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!(s.mean().unwrap().is_nan());
        assert!(s.variance().unwrap().is_nan());
    }

    #[test]
    fn min_max_empty_returns_none() {
        let s = RunningMoments::<2>::default();
        assert_eq!(Min::min(&s), None);
        assert_eq!(Max::max(&s), None);
    }

    #[test]
    fn min_max_populated() {
        let s = [3.0_f64, 1.0, 4.0, 1.0, 5.0]
            .iter().copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert_eq!(Min::min(&s), Some(1.0));
        assert_eq!(Max::max(&s), Some(5.0));
    }
}
