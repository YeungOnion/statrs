use crate::Sealed;
use crate::experimental_api::fold::Accumulate;
use crate::experimental_api::{Moments, PopulationMoments, Skewness};

pub struct RunningMoments<const ORDER: usize> {
    pub count: u64,
    m: [f64; ORDER],
    c: [f64; ORDER],
}

impl<const ORDER: usize> Default for RunningMoments<ORDER> {
    fn default() -> Self {
        Self {
            count: 0,
            m: [0.0; ORDER],
            c: [0.0; ORDER],
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

        // Welford / Pébaÿ (2008) central moment update.
        // m[0] = running mean (M1)
        // m[1] = sum of squared deviations (M2), if ORDER >= 2
        // m[2] = sum of cubed deviations (M3), if ORDER >= 3
        //
        // Update order: M3 before M2 before M1; each increment uses the
        // previous observation's lower-order accumulators.
        let delta = x - self.m[0];
        let delta_n = delta / n;
        let new_mean = self.m[0] + delta_n;
        let delta2 = x - new_mean;

        if let Some(&old_m2) = self.m.get(1) {
            if self.m.get(2).is_some() {
                let delta_n2 = delta_n * delta_n;
                self.compensated_add(2, delta * delta_n2 * (n - 1.0) * (n - 2.0) - 3.0 * delta_n * old_m2);
            }
            self.compensated_add(1, delta * delta2);
        }

        self.m[0] = new_mean;
        self
    }

    // Kahan-compensated add to m[k]: caller has verified k < ORDER.
    fn compensated_add(&mut self, k: usize, inc: f64) {
        let y = inc - self.c[k];
        let t = self.m[k] + y;
        self.c[k] = (t - self.m[k]) - y;
        self.m[k] = t;
    }

    pub fn from_iter<I: IntoIterator<Item = f64>>(iter: I) -> Self {
        iter.into_iter().fold(Self::default(), Self::push)
    }
}

impl Sealed for RunningMoments<2> {}
impl Sealed for RunningMoments<3> {}

impl Moments for RunningMoments<3> {
    fn mean(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.m[0])
        }
    }

    fn variance(&self) -> Option<f64> {
        if self.count < 2 {
            None
        } else {
            Some(self.m[1] / (self.count - 1) as f64)
        }
    }
}

impl Skewness for RunningMoments<3> {
    fn skewness(&self) -> Option<f64> {
        if self.count < 2 {
            return None;
        }
        let n = self.count as f64;
        let m2_mean = self.m[1] / n;
        let m3_mean = self.m[2] / n;
        let denom = m2_mean.powf(1.5);
        if denom == 0.0 {
            Some(0.0)
        } else {
            Some(m3_mean / denom)
        }
    }
}

impl Moments for RunningMoments<2> {
    fn mean(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.m[0])
        }
    }

    fn variance(&self) -> Option<f64> {
        if self.count < 2 {
            None
        } else {
            Some(self.m[1] / (self.count - 1) as f64)
        }
    }
}

impl PopulationMoments for RunningMoments<2> {
    fn population_variance(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.m[1] / self.count as f64)
        }
    }
}

impl PopulationMoments for RunningMoments<3> {
    fn population_variance(&self) -> Option<f64> {
        if self.count == 0 {
            None
        } else {
            Some(self.m[1] / self.count as f64)
        }
    }
}

impl<const ORDER: usize> Accumulate for RunningMoments<ORDER> {
    fn push(self, x: f64) -> Self {
        RunningMoments::push(self, x)
    }
}

/// Single-pass mean accumulator. Alias for [`RunningMoments<2>`] — variance is
/// computed alongside mean at no extra cost via Welford's algorithm.
pub type MeanAccum = RunningMoments<2>;

/// Single-pass N-variable sample covariance matrix accumulator.
///
/// `push` updates only the upper triangle of the cross-product matrix;
/// `finalize` mirrors it to the lower triangle and scales by `n − 1`.
///
/// Compose with [`RunningMoments`] via `fold` to get means and the covariance
/// matrix in one pass.
pub struct RunningCov<const N: usize> {
    pub count: u64,
    mean: [f64; N],
    /// Upper triangle of the running cross-product matrix (Welford online).
    /// `c[i][j]` is only maintained for `j >= i`; lower triangle is zeroed
    /// until `finalize` fills it in.
    c: [[f64; N]; N],
}

impl<const N: usize> Default for RunningCov<N> {
    fn default() -> Self {
        Self {
            count: 0,
            mean: [0.0; N],
            c: [[0.0; N]; N],
        }
    }
}

impl<const N: usize> RunningCov<N> {
    pub fn push(mut self, x: [f64; N]) -> Self {
        self.count += 1;
        let n = self.count as f64;

        let mut delta = [0.0; N];
        for i in 0..N {
            delta[i] = x[i] - self.mean[i];
            self.mean[i] += delta[i] / n;
        }

        // delta2[j] = x[j] - new_mean[j]; only upper triangle (j >= i)
        for i in 0..N {
            for j in i..N {
                self.c[i][j] += delta[i] * (x[j] - self.mean[j]);
            }
        }

        self
    }

    /// Returns the sample covariance matrix (normalised by `n − 1`), or
    /// `None` if fewer than two observations have been pushed.
    ///
    /// Mirrors the upper triangle to the lower triangle before returning.
    pub fn finalize(mut self) -> Option<[[f64; N]; N]> {
        if self.count < 2 {
            return None;
        }
        let denom = (self.count - 1) as f64;
        for i in 0..N {
            for j in i..N {
                self.c[i][j] /= denom;
                self.c[j][i] = self.c[i][j];
            }
        }
        Some(self.c)
    }
}

/// Type alias for a covariance accumulator over `N` variables.
pub type CovAccum<const N: usize> = RunningCov<N>;

/// Single-pass absolute minimum accumulator.
///
/// NaN poisons the state: any NaN observation causes [`abs_min`] to return `None`.
///
/// [`abs_min`]: AbsMinAccum::abs_min
#[derive(Default)]
pub struct AbsMinAccum(Option<f64>);

impl AbsMinAccum {
    pub fn abs_min(self) -> Option<f64> {
        self.0.filter(|v| !v.is_nan())
    }
}

impl Accumulate for AbsMinAccum {
    fn push(self, x: f64) -> Self {
        match self.0 {
            Some(v) if v.is_nan() => self,
            _ if x.is_nan() => Self(Some(f64::NAN)),
            Some(v) => Self(Some(v.min(x.abs()))),
            None => Self(Some(x.abs())),
        }
    }
}

/// Single-pass absolute maximum accumulator.
///
/// NaN poisons the state: any NaN observation causes [`abs_max`] to return `None`.
///
/// [`abs_max`]: AbsMaxAccum::abs_max
#[derive(Default)]
pub struct AbsMaxAccum(Option<f64>);

impl AbsMaxAccum {
    pub fn abs_max(self) -> Option<f64> {
        self.0.filter(|v| !v.is_nan())
    }
}

impl Accumulate for AbsMaxAccum {
    fn push(self, x: f64) -> Self {
        match self.0 {
            Some(v) if v.is_nan() => self,
            _ if x.is_nan() => Self(Some(f64::NAN)),
            Some(v) => Self(Some(v.max(x.abs()))),
            None => Self(Some(x.abs())),
        }
    }
}

/// Single-pass mean and variance accumulator.
pub type VarianceAccum = RunningMoments<2>;

/// Single-pass mean, variance, and skewness accumulator.
pub type SkewnessAccum = RunningMoments<3>;

#[cfg(test)]
mod tests {
    use super::*;

    type Accum = RunningMoments<2>;

    #[test]
    fn default_is_empty() {
        let s = Accum::default();
        assert_eq!(s.count, 0);
    }

    #[test]
    fn variance_and_min() {
        let reducer = |(s, m): (RunningMoments<_>, f64), x: f64| {
            if x.is_finite() {
                Ok((s.push(x), m.min(x)))
            } else {
                Err((x, x))
            }
        };
        let data = [1.0_f64, 2.0, 3.0].iter();
        let (var, m) = match data.copied().try_fold(
            (VarianceAccum::default(), f64::INFINITY),
            reducer as fn(_, _) -> _,
        ) {
            Ok((s, m)) => (s.variance().unwrap(), m),
            Err((x, y)) => (x, y),
        };
        assert_eq!(var, 1.0);
        assert_eq!(m, 1.0);
    }

    #[test]
    fn population_variance_empty_returns_none() {
        assert_eq!(RunningMoments::<2>::default().population_variance(), None);
    }

    #[test]
    fn population_variance_single_element() {
        // population variance of one element is 0 (no spread)
        let s = RunningMoments::<2>::default().push(5.0);
        assert_eq!(s.population_variance(), Some(0.0));
    }

    #[test]
    fn population_variance_known_dataset() {
        // [2,4,4,4,5,5,7,9]: population variance = M2/n = 32/8 = 4.0
        let data = [2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s = data
            .iter()
            .copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!((s.population_variance().unwrap() - 4.0).abs() < 1e-12);
        assert!((s.population_std_dev().unwrap() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn population_variance_less_than_sample_variance() {
        // population variance uses N, sample uses N-1 — population is always smaller
        let data = [1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let s = data
            .iter()
            .copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!(s.population_variance().unwrap() < s.variance().unwrap());
    }

    #[test]
    fn push_increments_count() {
        let s = Accum::default().push(3.0);
        assert_eq!(s.count, 1);
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
        let s = data
            .iter()
            .copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!((s.mean().unwrap() - 5.0).abs() < 1e-12);
        assert!((s.variance().unwrap() - 32.0 / 7.0).abs() < 1e-12);
        assert!((s.std_dev().unwrap() - (32.0_f64 / 7.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn moments_nan_propagates_as_some_nan() {
        let s = [1.0_f64, f64::NAN]
            .iter()
            .copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        assert!(s.mean().unwrap().is_nan());
        assert!(s.variance().unwrap().is_nan());
    }

    #[test]
    fn skewness_empty_returns_none() {
        assert_eq!(RunningMoments::<3>::default().skewness(), None);
    }

    #[test]
    fn skewness_single_returns_none() {
        assert_eq!(RunningMoments::<3>::default().push(1.0).skewness(), None);
    }

    #[test]
    fn skewness_symmetric_data_near_zero() {
        let data = [1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let s = data
            .iter()
            .copied()
            .fold(RunningMoments::<3>::default(), RunningMoments::push);
        assert!(s.skewness().unwrap().abs() < 1e-10);
    }

    #[test]
    fn skewness_known_dataset() {
        // [2,4,4,4,5,5,7,9]: population skewness = (M3/n) / (M2/n)^(3/2)
        // M2 = 32, M3 = 42, n = 8
        // skewness = (42/8) / (32/8)^(3/2) = 5.25 / 8.0 = 0.65625
        let data = [2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s = data
            .iter()
            .copied()
            .fold(RunningMoments::<3>::default(), RunningMoments::push);
        assert!((s.skewness().unwrap() - 0.65625).abs() < 1e-10);
    }

    #[test]
    fn mean_accum_and_variance_accum_are_same_type() {
        let _: fn(MeanAccum, f64) -> MeanAccum = VarianceAccum::push;
    }

    #[test]
    fn fold_ergonomics_example() {
        let data = [1.0_f64, 2.0, 3.0, 4.0, 5.0];
        let s = data
            .iter()
            .copied()
            .fold(VarianceAccum::default(), VarianceAccum::push);
        assert_eq!(s.count, 5);
        assert!((s.mean().unwrap() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn skewness_accum_alias_works() {
        let data = [1.0_f64, 2.0, 3.0];
        let s = data
            .iter()
            .copied()
            .fold(SkewnessAccum::default(), SkewnessAccum::push);
        assert_eq!(s.count, 3);
        assert!(s.skewness().is_some());
    }

    #[test]
    fn skewness_moments_consistent_with_order2() {
        let data = [2.0_f64, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0];
        let s2 = data
            .iter()
            .copied()
            .fold(RunningMoments::<2>::default(), RunningMoments::push);
        let s3 = data
            .iter()
            .copied()
            .fold(RunningMoments::<3>::default(), RunningMoments::push);
        assert!((s2.mean().unwrap() - s3.mean().unwrap()).abs() < 1e-12);
        assert!((s2.variance().unwrap() - s3.variance().unwrap()).abs() < 1e-12);
    }
}

#[cfg(test)]
mod cov_tests {
    use super::*;

    #[test]
    fn default_count_is_zero() {
        assert_eq!(RunningCov::<2>::default().count, 0);
    }

    #[test]
    fn push_increments_count() {
        assert_eq!(RunningCov::<2>::default().push([1.0, 2.0]).count, 1);
    }

    #[test]
    fn finalize_empty_is_none() {
        assert!(RunningCov::<2>::default().finalize().is_none());
    }

    #[test]
    fn finalize_single_obs_is_none() {
        assert!(
            RunningCov::<2>::default()
                .push([1.0, 2.0])
                .finalize()
                .is_none()
        );
    }

    #[test]
    fn finalize_n2_positive_correlation() {
        // [(1,2),(3,4),(5,6)]: means=[3,4], cov matrix = [[4,4],[4,4]]
        let s = [[1.0_f64, 2.0], [3.0, 4.0], [5.0, 6.0]]
            .into_iter()
            .fold(RunningCov::<2>::default(), RunningCov::push);
        let mat = s.finalize().unwrap();
        assert!((mat[0][0] - 4.0).abs() < 1e-12);
        assert!((mat[0][1] - 4.0).abs() < 1e-12);
        assert!((mat[1][0] - 4.0).abs() < 1e-12);
        assert!((mat[1][1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn finalize_n2_negative_correlation() {
        // [(1,6),(3,4),(5,2)]: means=[3,4], cov matrix = [[4,-4],[-4,4]]
        let s = [[1.0_f64, 6.0], [3.0, 4.0], [5.0, 2.0]]
            .into_iter()
            .fold(RunningCov::<2>::default(), RunningCov::push);
        let mat = s.finalize().unwrap();
        assert!((mat[0][0] - 4.0).abs() < 1e-12);
        assert!((mat[0][1] + 4.0).abs() < 1e-12);
        assert!((mat[1][0] + 4.0).abs() < 1e-12);
        assert!((mat[1][1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn finalize_is_symmetric() {
        let s = [[1.0_f64, 6.0], [3.0, 4.0], [5.0, 2.0]]
            .into_iter()
            .fold(RunningCov::<2>::default(), RunningCov::push);
        let mat = s.finalize().unwrap();
        assert_eq!(mat[0][1], mat[1][0]);
    }

    #[test]
    fn n1_covariance_matches_variance() {
        // N=1: cov matrix [[variance]]; [1,2,3] has sample variance 1.0
        let s = [[1.0_f64], [2.0], [3.0]]
            .into_iter()
            .fold(RunningCov::<1>::default(), RunningCov::push);
        let mat = s.finalize().unwrap();
        assert!((mat[0][0] - 1.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod accumulate_tests {
    use super::*;
    use crate::experimental_api::fold::Accumulate;

    #[test]
    fn running_moments_impl_accumulate() {
        let s: RunningMoments<2> = [1.0_f64, 2.0, 3.0]
            .iter()
            .copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!(s.mean(), Some(2.0));
    }

    #[test]
    fn abs_min_empty_is_none() {
        assert_eq!(AbsMinAccum::default().abs_min(), None);
    }

    #[test]
    fn abs_min_basic() {
        let acc = AbsMinAccum::default().push(3.0).push(-1.0).push(4.0);
        assert_eq!(acc.abs_min(), Some(1.0));
    }

    #[test]
    fn abs_min_nan_poisons() {
        let acc = AbsMinAccum::default().push(1.0).push(f64::NAN).push(0.5);
        assert_eq!(acc.abs_min(), None);
    }

    #[test]
    fn abs_max_empty_is_none() {
        assert_eq!(AbsMaxAccum::default().abs_max(), None);
    }

    #[test]
    fn abs_max_basic() {
        let acc = AbsMaxAccum::default().push(-5.0).push(2.0);
        assert_eq!(acc.abs_max(), Some(5.0));
    }

    #[test]
    fn abs_max_nan_poisons() {
        let acc = AbsMaxAccum::default().push(10.0).push(f64::NAN).push(20.0);
        assert_eq!(acc.abs_max(), None);
    }

    #[test]
    fn tuple_composition_variance_and_abs_min() {
        let data = [3.0_f64, -1.0, 4.0, -1.0, 5.0];
        let (var, mn): (VarianceAccum, AbsMinAccum) = data
            .iter()
            .copied()
            .fold(Default::default(), Accumulate::push);
        assert_eq!(mn.abs_min(), Some(1.0));
        assert!(var.variance().is_some());
    }
}
