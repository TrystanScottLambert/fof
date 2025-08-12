use std::collections::HashSet;
use rustc_hash::FxHashSet;

use crate::spherical_trig_funcs::convert_equitorial_to_cartesian;
use kiddo::{ImmutableKdTree, SquaredEuclidean};
use rayon::prelude::*;

// Combined function that finds indices and removes target in one go
fn find_indices_in_range(
    sorted_array: &[f64],
    argsort_array: &[usize],
    low_lim: f64,
    up_lim: f64,
    exclude_idx: usize,
) -> Vec<usize> {
    let start_idx = sorted_array.partition_point(|&x| x < low_lim);
    if start_idx >= sorted_array.len() {
        return Vec::new();
    }

    let end_idx = sorted_array.partition_point(|&x| x <= up_lim);
    if start_idx >= end_idx {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(end_idx - start_idx - 1);
    for &idx in &argsort_array[start_idx..end_idx] {
        if idx > exclude_idx {
            result.push(idx);
        }
    }
    result
}

fn argsort<T: PartialOrd>(data: &[T]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..data.len()).collect();
    idx.sort_by(|&i, &j| data[i].partial_cmp(&data[j]).unwrap());
    idx
}

/// Copy of the c++ implementation of the original gama-group finder that we can test against.
pub fn ffl1(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    // Convert (RA, Dec, dist) to 3D Cartesian coordinates.
    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();

    let mut ind = Vec::new();

    for i in 0..(n - 1) {
        for j in (i + 1)..n {
            let ztemp = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
            let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

            if zrad <= ztemp {
                let bgal2 = ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                let radproj = (0..3)
                    .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                    .sum::<f64>();

                if radproj <= bgal2 {
                    ind.push((i, j));
                }
            }
        }
    }
    ind
}

// Faster version where we presort the radial direction. Used for testing not production.
pub fn find_links_just_z(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {
    let n = ra_array.len();

    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();

    let max_los_ll = linking_lengths_los.iter().cloned().fold(f64::NAN, f64::max);

    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let dist_argsort = argsort(&comoving_distances);

    // Parallel outer loop
    let results: Vec<(usize, usize)> = (0..(n - 1))
        .into_par_iter()
        .map(|i| {
            let dist_i = comoving_distances[i];
            let lower_lim = dist_i - (max_los_ll + linking_lengths_los[i]) * 0.5;
            let upper_lim = dist_i + (max_los_ll + linking_lengths_los[i]) * 0.5;
            let possible_los_idx =
                find_indices_in_range(&sorted_distances, &dist_argsort, lower_lim, upper_lim, i);

            let mut local_pairs = Vec::new();

            for j in possible_los_idx {
                let average_los_ll = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
                let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

                if zrad <= average_los_ll {
                    let bgal2 = ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                    let radproj = (0..3)
                        .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                        .sum::<f64>();

                    if radproj <= bgal2 {
                        local_pairs.push((i, j));
                    }
                }
            }

            local_pairs
        })
        .flatten()
        .collect();

    results
}

/// Find all the connections between all galaxies in a redshift survey.
/// Returns a vector of tuples of (i, j) for galaxy i and j respectively.
/// This can then be used to construct the group catalog.
pub fn find_links(
    ra_array: Vec<f64>,
    dec_array: Vec<f64>,
    comoving_distances: Vec<f64>,
    linking_lengths_pos: Vec<f64>,
    linking_lengths_los: Vec<f64>,
) -> Vec<(usize, usize)> {

    let n = {
        let sizes = vec![ra_array.len(), dec_array.len(), comoving_distances.len(), linking_lengths_los.len(), linking_lengths_pos.len()]
            .into_iter()
            .collect::<HashSet<_>>();
        if sizes.len() != 1 {
            panic!("Inputs are not of the same size")
        }
        sizes.into_iter().next().unwrap()
    };

    let coords: Vec<[f64; 3]> = (0..n)
        .map(|i| convert_equitorial_to_cartesian(&ra_array[i], &dec_array[i]))
        .collect();

    let global_tree = ImmutableKdTree::new_from_slice(&coords);

    let max_los_ll = linking_lengths_los.iter().cloned().fold(f64::NAN, f64::max);

    let mut sorted_distances = comoving_distances.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let dist_argsort = argsort(&comoving_distances);

    // Parallel outer loop
    let results: Vec<(usize, usize)> = (0..(n - 1))
        .into_par_iter()
        // .with_min_len(1023)
        .map(|i| {
            let point = [coords[i][0], coords[i][1], coords[i][2]];
            let dist_i = comoving_distances[i];
            let lower_lim = dist_i - (max_los_ll + linking_lengths_los[i]) * 0.5;
            let upper_lim = dist_i + (max_los_ll + linking_lengths_los[i]) * 0.5;

            let possible_los_idx =
                find_indices_in_range(&sorted_distances, &dist_argsort, lower_lim, upper_lim, i);

            let max_local_pos_ll = {
                let max = possible_los_idx
                    .iter()
                    .map(|&idx| linking_lengths_pos[idx])
                    .reduce(f64::max);
                match max {
                    None => 0.,
                    Some(val) => (val + linking_lengths_pos[i]) * 0.5,
                }
            };

            let global_search =
                global_tree.within_unsorted::<SquaredEuclidean>(&point, max_local_pos_ll.powi(2));
            let global_set: FxHashSet<usize> = global_search
                .iter()
                .filter(|&n| n.item > (i as u64))
                .map(|n| n.item as usize)
                .collect();

            let mut local_pairs = Vec::new();

            possible_los_idx
                .into_iter()
                .filter(|j| global_set.contains(j))
                .for_each(|j| {
                    let average_los_ll = (linking_lengths_los[i] + linking_lengths_los[j]) * 0.5;
                    let zrad = (comoving_distances[i] - comoving_distances[j]).abs();

                    if zrad <= average_los_ll {
                        let bgal2 =
                            ((linking_lengths_pos[i] + linking_lengths_pos[j]) * 0.5).powi(2);

                        let radproj = (0..3)
                            .map(|k| (coords[i][k] - coords[j][k]).powi(2))
                            .sum::<f64>();

                        if radproj <= bgal2 {
                            local_pairs.push((i, j));
                        }
                    }
                });

            local_pairs
        })
        .flatten()
        .collect();

    results
}

#[cfg(test)]
mod tests {
    use std::iter::zip;

    use super::*;

    #[test]
    fn faster_links() {
        // These functions should be better benchmarked in the R code reading in the
        // old gama data.

        // Doesn't really matter what we put here. They have to match FFL1 exactly.
        // 3 2 1 group set up.
        let ra_array = vec![20., 20., 20., -20., -20., 0.];
        let dec_array = vec![-50., -50., -50., 50., 50., 75.];
        let comoving_distances = vec![100., 100., 100., 20., 20., 200.];
        let linking_lengths_pos = vec![0.06; 6];
        let linking_lengths_los = vec![18.; 6];

        let classic_links = ffl1(
            ra_array.clone(),
            dec_array.clone(),
            comoving_distances.clone(),
            linking_lengths_pos.clone(),
            linking_lengths_los.clone(),
        );
        let new_links = find_links(
            ra_array.clone(),
            dec_array.clone(),
            comoving_distances.clone(),
            linking_lengths_pos.clone(),
            linking_lengths_los.clone(),
        );
        for (new, classic) in zip(new_links, classic_links) {
            assert_eq!(new.0, classic.0);
            assert_eq!(new.1, classic.1);
        }
    }

    #[test]
    fn argsort_simple() {
        let data = [5., 3., 4., 1., 3.];
        let answer: [usize; 5] = [3, 1, 4, 2, 0];
        let result = argsort(&data);
        for (a, r) in zip(answer, result) {
            assert_eq!(a, r)
        }
    }

    #[test]
    fn test_find_indices_in_range() {
        let _data = [5., 3., 4., 1., 3., 0., 10.];
        let ans = [4];
        let data_sorted = [0., 1., 3., 3., 4., 5., 10.];
        let argsort = [5, 3, 1, 4, 2, 0, 6];
        let lower_lim = 3.;
        let upper_lim = 6.;

        let idx = find_indices_in_range(&data_sorted, &argsort, lower_lim, upper_lim, 2);
        assert_eq!(idx, &ans);
    }
}
