pub fn quadratic_mean(data: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut count = 0usize;
    for x in data {
        if x.is_nan() {
            return f64::NAN;
        }
        sum += x * x;
        count += 1;
    }
    if count == 0 {
        return f64::NAN;
    }
    (sum / count as f64).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_is_nan() {
        assert!(quadratic_mean([] as [f64; 0]).is_nan());
    }

    #[test]
    fn nan_entry_propagates() {
        assert!(quadratic_mean([1.0, f64::NAN, 2.0]).is_nan());
    }

    #[test]
    fn known_value() {
        let result = quadratic_mean([0.0, 3.0, -2.0]);
        assert!((result - 2.08167_f64).abs() < 1e-4);
    }
}
