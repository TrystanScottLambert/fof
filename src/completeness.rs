use kiddo::SquaredEuclidean;

use crate::spherical_trig_funcs::{chord_distance, convert_equitorial_to_cartesian};
use crate::tree::build_kd_tree;

#[derive(Clone, Debug)]
pub struct PositionCatalog {
    pub ra_deg: Vec<f64>,
    pub dec_deg: Vec<f64>
}

pub fn count_galaxies_around_point(
    catalog: PositionCatalog,
    ra_point_deg: f64,
    dec_point_deg: f64,
    angular_difference_deg: f64,
) -> usize {
    let chord_length = chord_distance(angular_difference_deg);
    let tree = build_kd_tree(catalog.ra_deg, catalog.dec_deg);
    let point_cartesian = convert_equitorial_to_cartesian(&ra_point_deg, &dec_point_deg);
    let within = tree.within::<SquaredEuclidean>(&point_cartesian, chord_length.powi(2));
    within.len()
}

pub fn count_galaxies_around_points(
    catalog: PositionCatalog,
    evaluate_catalog: PositionCatalog,
    angular_differences_deg: Vec<f64>,
) -> Vec<usize> {
    let chord_lengths = angular_differences_deg
        .iter()
        .map(|&ang| chord_distance(ang));

    let tree = build_kd_tree(catalog.ra_deg, catalog.dec_deg);
    let points_cartesian: Vec<[f64; 3]> = evaluate_catalog.ra_deg
        .iter()
        .zip(evaluate_catalog.dec_deg)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(ra, &dec))
        .collect();

    points_cartesian
        .iter()
        .zip(chord_lengths)
        .map(|(point, chord)| tree.within::<SquaredEuclidean>(point, chord.powi(2)).len())
        .collect()
}

pub fn calculate_completeness(observed_catalog: PositionCatalog, target_catalog: PositionCatalog, evaluate_catalog: PositionCatalog, angular_radius_deg: Vec<f64>) -> Vec<f64> {
    let number_observed = count_galaxies_around_points(observed_catalog, evaluate_catalog.clone(), angular_radius_deg.clone());
    let number_target = count_galaxies_around_points(target_catalog, evaluate_catalog, angular_radius_deg);
    number_observed.iter()
        .zip(number_target.iter())
        .map(|(obs, tar)| (*obs as f64)/(*tar as f64))
        .collect::<Vec<f64>>()
}


#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn test_counting_around_a_single_point() {
        let ras = vec![120., 120.9999, 120., 180.];
        let decs = vec![0., 0., 0., -45.];
        let cat = PositionCatalog {ra_deg: ras, dec_deg: decs};
        let result_3 = count_galaxies_around_point(cat.clone(), 120., 0., 1.);
        let result_2 = count_galaxies_around_point(cat.clone(), 120., 0., 0.1);
        let result_1 = count_galaxies_around_point(cat.clone(), 179., -46., 5.);

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
        let catalog = PositionCatalog {ra_deg: ras, dec_deg: decs};
        let eval_catalog = PositionCatalog {ra_deg: ra_points, dec_deg: dec_points};
        let ang_seps = vec![1., 0.1, 5.];
        let results = count_galaxies_around_points(catalog, eval_catalog, ang_seps);
        let answers = [3, 2, 1];
        for (res, ans) in zip(results, answers) {
            assert_eq!(res, ans)
        }
    }

    #[test]
    fn test_completeness() {
        let target = PositionCatalog {ra_deg: vec![20.;4], dec_deg: vec![-20.;4]};
        let observed = PositionCatalog {ra_deg: vec![20.;3], dec_deg: vec![-20.;3]};
        let eval = observed.clone();
        let ang_dist = vec![0.1;3];

        let result = calculate_completeness(observed, target, eval, ang_dist);
        for (res, ans) in zip(result, vec![0.75; 3]) {
            assert_eq!(res, ans)
        }

    }

}
