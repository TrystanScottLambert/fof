use crate::spherical_trig_funcs::{chord_distance, convert_equitorial_to_cartesian};
use kiddo::ImmutableKdTree;
use kiddo::SquaredEuclidean;

pub struct Point {
    pub ra_deg: f64,
    pub dec_deg: f64,
}

pub fn build_kd_tree(ra_array_deg: Vec<f64>, dec_array_deg: Vec<f64>) -> ImmutableKdTree<f64, 3> {
    let entries: Vec<[f64; 3]> = ra_array_deg
        .iter()
        .zip(dec_array_deg)
        .map(|(ra, dec)| convert_equitorial_to_cartesian(ra, &dec))
        .collect();

    ImmutableKdTree::new_from_slice(&entries)
}

pub fn find_idx_within(tree: ImmutableKdTree<f64, 3>, point: Point, angular_difference_deg: f64) -> Vec<u64> {
    let point_cartesian = convert_equitorial_to_cartesian(&point.ra_deg, &point.dec_deg);
    let dist = chord_distance(angular_difference_deg).powi(2); // squared chord distance.
    let neighbors = tree.within::<SquaredEuclidean>(&point_cartesian, dist);
    neighbors.iter().map(|nn| nn.item).collect()
}

#[cfg(test)]
mod tests {

    use std::iter::zip;

    use super::*;

    #[test]
    fn test_kd_tree_building() {
        let ras = vec![0., 23., 180.];
        let decs = vec![0., -10., 20.];
        let tree = build_kd_tree(ras.clone(), decs.clone());

        let test_point = convert_equitorial_to_cartesian(&ras[1], &decs[1]);
        let result = tree.within::<SquaredEuclidean>(&test_point, 1e-10);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_finding_within() {
        let ras = vec![120., 120.9999, 120., 180.];
        let decs = vec![0., 0., 0., -45.];
        let tree = build_kd_tree(ras, decs);
        let point = Point { ra_deg: 120., dec_deg: 0. };
        let result = find_idx_within(tree, point, 1.);
        assert_eq!(result.len(), 3);
        let answers = [0, 2, 1];
        for (res, ans) in zip(result, answers) {
            assert_eq!(res, ans)
        }
    }

}