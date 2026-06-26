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
}
