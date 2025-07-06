/// Quantile caclulation. Produces the quantile of the given probability.
/// This is a clone of the quantile function in R.
///
pub fn quantile_interpolated(sorted: &[f64], quantile: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return f64::NAN;
    }
    let h = (n - 1) as f64 * quantile;
    let i = h.floor() as usize;
    let frac = h - i as f64;

    if i + 1 < n {
        sorted[i] * (1.0 - frac) + sorted[i + 1] * frac
    } else {
        sorted[i]
    }
}

/// Harmonic mean calculation. The reciprocal of the arithmetic mean of the reciprocals.
/// Meant to replicate the pscyh::harmonic.mean in R.
///
/// #Examples
/// ```
/// use fof::stats::harmonic_mean;
/// let values = vec![1., 4., 4.];
/// let answer = 2.;
/// let result = harmonic_mean(values);
/// assert_eq!(result, answer);
/// ```
pub fn harmonic_mean(values: Vec<f64>) -> f64 {
    let n = values.iter().len() as f64;
    let summation = values.iter().map(|v| 1. / v).sum::<f64>();
    n / summation
}

/// Compute the arithmetic mean.
pub fn mean(values: Vec<f64>) -> f64 {
    let n = values.iter().len() as f64;
    let summation = values.iter().sum::<f64>();
    summation / n
}

/// Compute the sample median.
pub fn median(mut values: Vec<f64>) -> f64 {
    let len = values.len();
    if len == 0 {
        panic!("Cannot compute median of empty vector");
    }

    values.sort_by(|a, b| a.partial_cmp(b).unwrap());

    if len % 2 == 1 {
        values[len / 2]
    } else {
        let mid = len / 2;
        (values[mid - 1] + values[mid]) / 2.0
    }
}

/// Direct copy of the density function in R.
pub fn density(data: &[f64], binwidth: f64, from: f64, to: f64, n: usize) -> (Vec<f64>, Vec<f64>) {
    let step = (to - from) / (n as f64);
    let mut x_vals = Vec::with_capacity(n);
    let mut y_vals = Vec::with_capacity(n);

    for i in 0..n {
        let x = from + i as f64 * step;
        x_vals.push(x);

        let count = data
            .iter()
            .filter(|&&d| (x - binwidth / 2.0..x + binwidth / 2.0).contains(&d))
            .count();

        let density = count as f64 / (data.len() as f64 * binwidth);
        y_vals.push(density)
    }
    (x_vals, y_vals)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn comparing_quantile_to_r() {
        let sorted_array = [1., 2., 3., 4., 5.];
        assert!((quantile_interpolated(&sorted_array, 0.3) - 2.2).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 0.5) - 3.).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 0.77) - 4.08).abs() < 1e-6);
        assert!((quantile_interpolated(&sorted_array, 1.) - 5.).abs() < 1e-6);
    }

    #[test]
    fn comparing_to_psych_in_r_harmonic_mean() {
        // Checking that arrays with a zero are zero.
        let values = vec![0., 1., 2., 3.];
        let values_float = vec![0., 0.1, 0., 0.2, 0.3, 0.5];
        let answer = 0.;
        let result = harmonic_mean(values);
        let result_float = harmonic_mean(values_float);
        assert_eq!(result, answer);
        assert_eq!(result_float, answer);

        // Example often used.
        let values = vec![1., 4., 4.];
        let answer = 2.;
        let result = harmonic_mean(values);
        assert_eq!(result, answer);

        // Testing on more realistic floats.
        let values = vec![0.2, 0.3, 0.21, 0.22, 0.31, 0.5, 0.1, 0.12];
        let answers = 0.1941755;
        let result = harmonic_mean(values);
        assert!((answers - result).abs() < 1e-5)
    }
    #[test]
    fn testing_the_mean_function() {
        let values_1 = vec![0., 1., 2., 3.];
        let values_2 = vec![0., 0.1, 0.2, 0.3, 0.4];
        let values_3 = vec![-0.2, -0.23, -0.5, 0.1, 0.2];

        let answers_1 = 1.5;
        let answers_2 = 0.2;
        let answers_3 = -0.126;

        let results_1 = mean(values_1);
        let results_2 = mean(values_2);
        let results_3 = mean(values_3);

        assert!((results_1 - answers_1).abs() < 1e-5);
        assert!((results_2 - answers_2).abs() < 1e-5);
        assert!((results_3 - answers_3).abs() < 1e-5);
    }

    #[test]
    fn testing_the_median_function() {
        let values_1 = vec![0., 1., 2., 3.];
        let values_2 = vec![0., 0.1, 0.2, 0.3, 0.4];
        let values_3 = vec![-0.2, -0.23, -0.5, 0.1, 0.2];
        let values_4 = vec![0.1, 0.2, 0.1, 0.7, 0.2, 0.1, 0.8, 0.6, -0.1];

        let answers_1 = 1.5;
        let answers_2 = 0.2;
        let answers_3 = -0.2;
        let answers_4 = 0.2;

        let results_1 = median(values_1);
        let results_2 = median(values_2);
        let results_3 = median(values_3);
        let results_4 = median(values_4);

        assert!((results_1 - answers_1).abs() < 1e-5);
        assert!((results_2 - answers_2).abs() < 1e-5);
        assert!((results_3 - answers_3).abs() < 1e-5);
        assert!((results_4 - answers_4).abs() < 1e-5);
    }
}
