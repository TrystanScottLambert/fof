use kiddo::ImmutableKdTree;
use kiddo::SquaredEuclidean;

use crate::spherical_trig_funcs::{chord_distance, convert_equitorial_to_cartesian};

fn build_kd_tree(ra_array_deg: Vec<f64>, dec_array_deg: Vec<f64>) -> ImmutableKdTree<f64, 3> {
    let entries: Vec<[f64; 3]> = ra_array_deg
        .iter()
        .zip(dec_array_deg)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(ra, &dec))
        .collect();

    ImmutableKdTree::new_from_slice(&entries)
}

fn count_galaxies_around_point(
    ra_array_deg: Vec<f64>,
    dec_array_deg: Vec<f64>,
    ra_point_deg: f64,
    dec_point_deg: f64,
    angular_difference_deg: f64,
) -> usize {
    let chord_length = chord_distance(angular_difference_deg);
    let tree = build_kd_tree(ra_array_deg, dec_array_deg);
    let point_cartesian = convert_equitorial_to_cartesian(&ra_point_deg, &dec_point_deg);
    let within = tree.within::<SquaredEuclidean>(
        &point_cartesian,
        chord_length.powi(2),
    );
    within.len()
}

fn count_galaxies_around_points(
    ra_array_deg: Vec<f64>,
    dec_array_deg: Vec<f64>,
    ra_points_deg: Vec<f64>,
    dec_points_deg: Vec<f64>,
    angular_differences_deg: Vec<f64>,
) -> Vec<usize> {
    let chord_lengths = angular_differences_deg.iter().map(|&ang| chord_distance(ang));

    let tree = build_kd_tree(ra_array_deg, dec_array_deg);
    let points_cartesian: Vec<[f64; 3]> = ra_points_deg
        .iter()
        .zip(dec_points_deg)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(ra, &dec))
        .collect();

    points_cartesian
        .iter()
        .zip(chord_lengths)
        .map(|(point, chord)| {
            tree.within::<SquaredEuclidean>(
                point,
                chord.powi(2),
            )
            .len()
        })
        .collect()
    
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_kd_tree_building() {
        let ras = vec![0., 23., 180.];
        let decs  = vec![0., -10., 20.];
        let tree = build_kd_tree(ras, decs);
    }

    #[test]
    fn test_counting_around_a_single_point() {
        let ras = vec![120., 120.9999, 120., 180.];
        let decs = vec![0., 0., 0., -45.];
        let result_3 = count_galaxies_around_point(ras.clone(), decs.clone(), 120., 0., 1.);
        let result_2 = count_galaxies_around_point(ras.clone(), decs.clone(), 120., 0., 0.1);
        let result_1 = count_galaxies_around_point(ras, decs, 179., -46., 5.);

        assert_eq!(result_3, 3);
        assert_eq!(result_2, 2);
        assert_eq!(result_1, 1);

    }

    #[test]
    fn test_counting_around_multiple_points() {
        let ras = vec![120., 120.9999, 120., 180.];
        let decs = vec![0., 0., 0., -45.];
        let ra_points = vec![120., 120., 179.];
        let dec_points = vec![0., 0., -46.];
        let ang_seps = vec![1., 0.1, 5.];
        let results = count_galaxies_around_points(ras, decs, ra_points, dec_points, ang_seps);
        let answers = [3, 2, 1];
        for (res, ans) in zip(results, answers) {
            assert_eq!(res, ans)
        }

    }
}