use std::{f64::consts::PI, iter::zip};
use rayon::prelude::*;

use crate::constants::{G_MSOL_MPC_KMS2, SOLAR_MAG, SPEED_OF_LIGHT};
use crate::cosmology_funcs::Cosmology;
use crate::spherical_trig_funcs::{
    angular_separation, convert_cartesian_to_equitorial, convert_equitorial_to_cartesian,
    convert_equitorial_to_cartesian_scaled, euclidean_distance_3d,
};
use crate::stats::{mean, median, quantile_interpolated};

/// Calculating the total mass of the group from the R_g and 1d velocity dispersion
///
/// This is from Equation 8 of Tempel+2014 and assumes the viral theorem.
/// The gravitational_radius must be in Mpc and the los_velocity_dispersion is in km/s
/// Returns the mass in solar masses
fn calculate_total_mass(gravitational_radius: &f64, los_velocity_dispersion: &f64) -> f64 {
    2.325e12
        * gravitational_radius
        * ((3_f64.powf(1. / 3.)) * los_velocity_dispersion / 100.).powi(2)
}

/// A Group struct which stores the necessary values required for the group catalog.
pub struct Group {
    pub ra_members: Vec<f64>,
    pub dec_members: Vec<f64>,
    pub redshift_members: Vec<f64>,
    pub absolute_magnitude_members: Vec<f64>,
    pub velocity_errors: Vec<f64>,
}

impl Group {
    /// Determines the number of galaxies found in the group. Known as the multiplicity.
    pub fn multiplicity(&self) -> u32 {
        self.ra_members.len() as u32
    }

    /// The median of the redshfit of the members.
    pub fn median_redshift(&self) -> f64 {
        median(self.redshift_members.clone())
    }

    /// The comoving distance at the median redshift.
    pub fn median_distance(&self, cosmo: &Cosmology) -> f64 {
        cosmo.comoving_distance(self.median_redshift())
    }

    /// The velocity dispersion as calculated through the gapper method.
    /// The Sigma error would be in km/s.
    pub fn velocity_dispersion_gapper(&self) -> (f64, f64) {
        let sigma_err_squared = mean(self.velocity_errors.clone());
        let median_redshift = median(self.redshift_members.clone());
        let n = self.redshift_members.len();
        let nf64 = n as f64;

        let mut velocities: Vec<f64> = self
            .redshift_members
            .par_iter()
            .map(|z| (z * SPEED_OF_LIGHT) / (1. + median_redshift))
            .collect();

        velocities.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let gaps: Vec<f64> = velocities.windows(2).map(|w| w[1] - w[0]).collect();
        let i: Vec<usize> = (1..n).collect();
        let weights: Vec<usize> = i.par_iter().map(|val| val * (n - val)).collect();

        let sum: f64 = zip(weights, gaps).map(|(a, b)| a as f64 * b).sum();
        let sigma_gap: f64 = ((PI.sqrt()) / (nf64 * (nf64 - 1.))) * sum;
        let raw_dispersion_squared = (n as f64 * sigma_gap.powi(2)) / (nf64 - 1.);
        let dispersion = if raw_dispersion_squared > sigma_err_squared {
            (raw_dispersion_squared - sigma_err_squared).sqrt()
        } else {
            0.0
        };
        (dispersion, sigma_err_squared.sqrt())
    }

    /// Returns the index (in the original input arrays) of the final remaining object
    /// after iteratively removing the furthest objects from the flux-weighted center.
    pub fn calculate_iterative_center_idx(&self) -> usize {
        let coords_cartesian: Vec<[f64; 3]> =
            zip(self.ra_members.clone(), self.dec_members.clone())
                .map(|(ra, dec)| convert_equitorial_to_cartesian(&ra, &dec))
                .collect();

        let flux: Vec<f64> = self
            .absolute_magnitude_members
            .par_iter()
            .map(|mag| 10_f64.powf(-0.4 * mag))
            .collect();

        // Track original indices alongside coords and flux
        let mut temp_coords = coords_cartesian.clone();
        let mut temp_flux = flux.clone();
        let mut temp_indices: Vec<usize> = (0..coords_cartesian.len()).collect();

        while temp_flux.len() > 2 {
            let flux_sum: f64 = temp_flux.par_iter().cloned().sum();

            let center: [f64; 3] = (0..3)
                .map(|i| {
                    temp_coords
                        .iter()
                        .zip(temp_flux.iter())
                        .map(|(coord, &f)| coord[i] * f)
                        .sum::<f64>()
                        / flux_sum
                })
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap();

            let distances: Vec<f64> = temp_coords
                .par_iter()
                .map(|coord| euclidean_distance_3d(coord, &center))
                .collect();

            if let Some((max_idx, _)) = distances
                .par_iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            {
                temp_coords.remove(max_idx);
                temp_flux.remove(max_idx);
                temp_indices.remove(max_idx);
            } else {
                break;
            }
        }

        // Find the index of the remaining object with highest flux
        let max_flux_idx = temp_flux
            .par_iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .map(|(i, _)| i)
            .unwrap();

        // Return the original index
        temp_indices[max_flux_idx]
    }

    /// Dispersion of the plane of sky
    ///
    /// Calculating the dispersion of the plane of sky using Equation 4 from Tempel+2014
    /// Using the iterative center as the definition of center.
    /// In units of Mpc
    pub fn calculate_sky_distribution(&self, cosmo: &Cosmology) -> f64 {
        let iterative_idx = self.calculate_iterative_center_idx();
        let group_ra = self.ra_members[iterative_idx];
        let group_dec = self.dec_members[iterative_idx];
        let projected_distances: Vec<f64> = self
            .ra_members
            .par_iter()
            .zip(self.dec_members.par_iter())
            .map(|(ra, dec)| {
                angular_separation(ra, dec, &group_ra, &group_dec)
                    * 3600. // deg to arcseconds
                    * cosmo.kpc_per_arcsecond_comoving(self.median_redshift())
            })
            .collect();
        ((1. / ((self.ra_members.len() as f64 * 2.) * (1. + self.median_redshift()).powi(2)))
            * projected_distances
                .par_iter()
                .map(|kpc| (kpc / 1000.).powi(2)) // to Mpc
                .sum::<f64>())
        .sqrt()
    }

    /// Several radii measurements of the group.
    /// R50, Rsigma, R100. These are the radii that contain 50%, 68% and 100% of the galaxies.
    pub fn calculate_radius(
        &self,
        group_center_ra: f64,
        group_center_dec: f64,
        group_center_z: f64,
        cosmo: &Cosmology,
    ) -> [f64; 3] {
        let group_center_dist = cosmo.comoving_distance(group_center_z);
        let center = convert_equitorial_to_cartesian_scaled(
            group_center_ra,
            group_center_dec,
            group_center_dist,
        );

        let mut distances: Vec<f64> = self
            .ra_members
            .par_iter()
            .zip(&self.dec_members)
            .map(|(&ra, &dec)| {
                let pos = convert_equitorial_to_cartesian_scaled(ra, dec, group_center_dist);
                euclidean_distance_3d(&pos, &center)
            })
            .collect();

        distances.sort_by(|a, b| a.partial_cmp(b).unwrap());

        let q50 = quantile_interpolated(&distances, 0.5);
        let q68 = quantile_interpolated(&distances, 0.68);
        let q100 = *distances.last().unwrap();

        [q50, q68, q100]
    }

    pub fn calculate_center_of_light(&self) -> (f64, f64) {
        let fluxes: Vec<f64> = self
            .absolute_magnitude_members
            .par_iter()
            .map(|mag| 10.0_f64.powf(-0.4 * mag))
            .collect();
        let sum_flux: f64 = fluxes.par_iter().sum();
        let coords_cartesian: Vec<[f64; 3]> =
            zip(self.ra_members.clone(), self.dec_members.clone())
                .map(|(ra, dec)| convert_equitorial_to_cartesian(&ra, &dec))
                .collect();

        let weighted_x = coords_cartesian
            .par_iter()
            .zip(fluxes.par_iter())
            .map(|(coord, flux)| coord[0] * flux)
            .sum::<f64>()
            / sum_flux;

        let weighted_y = coords_cartesian
            .par_iter()
            .zip(fluxes.par_iter())
            .map(|(coord, flux)| coord[1] * flux)
            .sum::<f64>()
            / sum_flux;

        let weighted_z = coords_cartesian
            .par_iter()
            .zip(fluxes.par_iter())
            .map(|(coord, flux)| coord[2] * flux)
            .sum::<f64>()
            / sum_flux;

        let center = convert_cartesian_to_equitorial(&weighted_x, &weighted_y, &weighted_z);
        (center[0], center[1])
    }

    pub fn calculate_flux_weighted_redshift(&self) -> f64 {
        let fluxes: Vec<f64> = self
            .absolute_magnitude_members
            .par_iter()
            .map(|mag| 10_f64.powf(-0.4 * mag))
            .collect();
        let sum_flux: f64 = fluxes.par_iter().sum();

        self.redshift_members
            .par_iter()
            .zip(fluxes.par_iter())
            .map(|(red, flux)| red * flux)
            .sum::<f64>()
            / sum_flux
    }

    pub fn total_flux(&self) -> f64 {
        self.absolute_magnitude_members
            .par_iter()
            .map(|mag| 10.0_f64.powf(-0.4 * mag))
            .sum()
    }

    pub fn flux_proxy(&self) -> f64 {
        // Converting to Solar luminosities
        self.total_flux() * 10.0_f64.powf(0.4 * SOLAR_MAG)
    }

    pub fn total_absolute_magnitude(&self) -> f64 {
        -2.5 * self.total_flux().log10()
    }
}

/// Struct of the galaxy group catalog.
pub struct GroupedGalaxyCatalog {
    pub ra: Vec<f64>,
    pub dec: Vec<f64>,
    pub redshift: Vec<f64>,
    pub absolute_magnitudes: Vec<f64>,
    pub velocity_errors: Vec<f64>,
    pub group_ids: Vec<i32>,
}

impl GroupedGalaxyCatalog {
    /// Get all the group ids.
    fn get_unique_ids(&self) -> Vec<i32> {
        let mut ids = self.group_ids.clone();
        ids.sort();
        ids.dedup();
        ids
    }

    /// Calculates all the group properties.
    pub fn calculate_group_properties(&self, cosmo: &Cosmology) -> GroupCatalog {
        let unique_group_ids = self.get_unique_ids();
        let mut id_groups: Vec<i32> = Vec::new();
        let mut iterative_ras: Vec<f64> = Vec::new();
        let mut iterative_decs: Vec<f64> = Vec::new();
        let mut iterative_redshifts: Vec<f64> = Vec::new();
        let mut iterative_idxs: Vec<usize> = Vec::new();
        let mut median_redshifts: Vec<f64> = Vec::new();
        let mut distance_groups: Vec<f64> = Vec::new();
        let mut r50_groups: Vec<f64> = Vec::new();
        let mut r100_groups: Vec<f64> = Vec::new();
        let mut rsigma_groups: Vec<f64> = Vec::new();
        let mut multiplicity_groups: Vec<u32> = Vec::new();
        let mut velocity_dispersions: Vec<f64> = Vec::new();
        let mut velocity_dispersion_errs: Vec<f64> = Vec::new();
        let mut raw_masses: Vec<f64> = Vec::new();
        let mut estimated_masses: Vec<f64> = Vec::new();
        let mut bcg_idxs: Vec<usize> = Vec::new();
        let mut bcg_ras: Vec<f64> = Vec::new();
        let mut bcg_decs: Vec<f64> = Vec::new();
        let mut bcg_redshifts: Vec<f64> = Vec::new();
        let mut col_ras: Vec<f64> = Vec::new();
        let mut col_decs: Vec<f64> = Vec::new();
        let mut total_absolute_mags: Vec<f64> = Vec::new();
        let mut total_flux_proxies: Vec<f64> = Vec::new();

        for id in unique_group_ids {
            if id >= 0 {
                // ignoring singletons = -1
                let local_group_ids: Vec<i32> = self
                    .group_ids
                    .clone()
                    .into_iter()
                    .enumerate()
                    .filter_map(|(idx, i)| if i == id { Some(idx as i32) } else { None })
                    .collect();

                let local_ra: Vec<f64> = local_group_ids
                    .clone()
                    .into_iter()
                    .map(|i| *self.ra.get(i as usize).unwrap())
                    .collect();
                let local_dec: Vec<f64> = local_group_ids
                    .clone()
                    .into_iter()
                    .map(|i| *self.dec.get(i as usize).unwrap())
                    .collect();
                let local_z: Vec<f64> = local_group_ids
                    .clone()
                    .into_iter()
                    .map(|i| *self.redshift.get(i as usize).unwrap())
                    .collect();
                let local_mag: Vec<f64> = local_group_ids
                    .clone()
                    .into_iter()
                    .map(|i| *self.absolute_magnitudes.get(i as usize).unwrap())
                    .collect();
                let local_vel_errors: Vec<f64> = local_group_ids
                    .clone()
                    .into_iter()
                    .map(|i| *self.velocity_errors.get(i as usize).unwrap())
                    .collect();

                let local_bcg_id = local_mag
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .map(|(i, _)| i)
                    .unwrap();
                let global_bcg_idx = local_group_ids[local_bcg_id] as usize;
                let bcg_ra = local_ra[local_bcg_id];
                let bcg_dec = local_dec[local_bcg_id];
                let bcg_z = local_z[local_bcg_id];

                let local_group = Group {
                    ra_members: local_ra.clone(),
                    dec_members: local_dec.clone(),
                    redshift_members: local_z.clone(),
                    absolute_magnitude_members: local_mag,
                    velocity_errors: local_vel_errors,
                };

                let (velocity_disp, velocity_disp_err) = local_group.velocity_dispersion_gapper();

                //let (ra_group, dec_group) = local_group.calculate_iterative_center();
                let iterative_idx = local_group.calculate_iterative_center_idx();
                let (ra_group, dec_group, iterative_redshift) = (
                    local_ra[iterative_idx],
                    local_dec[iterative_idx],
                    local_z[iterative_idx],
                );
                let global_iterative_idx = local_group_ids[iterative_idx] as usize;

                let z_group = local_group.median_redshift();
                let [r50_group, rsimga_group, r100_group] =
                    local_group.calculate_radius(ra_group, dec_group, z_group, cosmo);

                let raw_mass = (r50_group * velocity_disp.powi(2)) / G_MSOL_MPC_KMS2;

                let (col_ra, col_dec) = local_group.calculate_center_of_light();

                let sky_disp = local_group.calculate_sky_distribution(cosmo);
                let grav_rad = 4.582 * sky_disp;
                let m_total = calculate_total_mass(&grav_rad, &velocity_disp);

                id_groups.push(id);
                iterative_ras.push(ra_group);
                iterative_decs.push(dec_group);
                iterative_redshifts.push(iterative_redshift);
                iterative_idxs.push(global_iterative_idx);
                median_redshifts.push(z_group);
                r50_groups.push(r50_group);
                r100_groups.push(r100_group);
                rsigma_groups.push(rsimga_group);
                distance_groups.push(local_group.median_distance(cosmo));
                multiplicity_groups.push(local_group.multiplicity());
                velocity_dispersions.push(velocity_disp);
                velocity_dispersion_errs.push(velocity_disp_err);
                raw_masses.push(raw_mass);
                estimated_masses.push(m_total);
                bcg_idxs.push(global_bcg_idx);
                bcg_ras.push(bcg_ra);
                bcg_decs.push(bcg_dec);
                bcg_redshifts.push(bcg_z);
                col_ras.push(col_ra);
                col_decs.push(col_dec);
                total_absolute_mags.push(local_group.total_absolute_magnitude());
                total_flux_proxies.push(local_group.flux_proxy());
            }
        }

        GroupCatalog {
            ids: id_groups,
            iter_ras: iterative_ras,
            iter_decs: iterative_decs,
            iter_redshifts: iterative_redshifts,
            iter_idxs: iterative_idxs,
            median_redshifts,
            distances: distance_groups,
            r50s: r50_groups,
            r100s: r100_groups,
            rsigmas: rsigma_groups,
            multiplicity: multiplicity_groups,
            velocity_dispersion_gap: velocity_dispersions,
            velocity_dispersion_gap_err: velocity_dispersion_errs,
            raw_masses,
            estimated_masses,
            bcg_idxs,
            bcg_ras,
            bcg_decs,
            bcg_redshifts,
            col_ras,
            col_decs,
            total_flux_proxies,
            total_absolute_mags,
        }
    }

    pub fn calculate_pair_properties(&self) -> PairCatalog {
        let unique_group_ids = self.get_unique_ids();
        let mut ids: Vec<i32> = Vec::new();
        let mut idx_1: Vec<i32> = Vec::new();
        let mut idx_2: Vec<i32> = Vec::new();
        let mut projected_separation: Vec<f64> = Vec::new();
        let mut velocity_separation: Vec<f64> = Vec::new();
        let mut ra_bar: Vec<f64> = Vec::new();
        let mut dec_bar: Vec<f64> = Vec::new();
        let mut redshift_bar: Vec<f64> = Vec::new();
        let mut total_absolute_mags: Vec<f64> = Vec::new();

        for id in unique_group_ids {
            if id >= 0 {
                // ignoring singletons = -1
                let local_group_ids: Vec<i32> = self
                    .group_ids
                    .clone()
                    .into_iter()
                    .enumerate()
                    .filter_map(|(idx, i)| if i == id { Some(idx as i32) } else { None })
                    .collect();

                if local_group_ids.len() == 2 {
                    let local_ra: Vec<f64> = local_group_ids
                        .clone()
                        .into_iter()
                        .map(|i| *self.ra.get(i as usize).unwrap())
                        .collect();
                    let local_dec: Vec<f64> = local_group_ids
                        .clone()
                        .into_iter()
                        .map(|i| *self.dec.get(i as usize).unwrap())
                        .collect();
                    let local_z: Vec<f64> = local_group_ids
                        .clone()
                        .into_iter()
                        .map(|i| *self.redshift.get(i as usize).unwrap())
                        .collect();
                    let local_mag: Vec<f64> = local_group_ids
                        .clone()
                        .into_iter()
                        .map(|i| *self.absolute_magnitudes.get(i as usize).unwrap())
                        .collect();

                    let local_group = Group {
                        ra_members: local_ra.clone(),
                        dec_members: local_dec.clone(),
                        redshift_members: local_z.clone(),
                        absolute_magnitude_members: local_mag,
                        velocity_errors: vec![50.; 1], // dummy variable
                    };

                    let (ra, dec) = local_group.calculate_center_of_light();

                    ids.push(id);
                    idx_1.push(local_group_ids[0]);
                    idx_2.push(local_group_ids[1]);
                    projected_separation.push(angular_separation(
                        &local_ra[0],
                        &local_dec[0],
                        &local_ra[1],
                        &local_dec[1],
                    ));
                    velocity_separation.push((local_z[0] - local_z[1]).abs());
                    ra_bar.push(ra);
                    dec_bar.push(dec);
                    redshift_bar.push(local_group.calculate_flux_weighted_redshift());
                    total_absolute_mags.push(local_group.total_absolute_magnitude())
                }
            }
        }
        PairCatalog {
            ids,
            idx_1,
            idx_2,
            projected_separation,
            velocity_separation,
            ra_bar,
            dec_bar,
            redshift_bar,
            total_absolute_mags,
        }
    }
}

/// Struct which represents the group catalog.
pub struct GroupCatalog {
    pub ids: Vec<i32>,
    pub iter_ras: Vec<f64>,
    pub iter_decs: Vec<f64>,
    pub iter_redshifts: Vec<f64>,
    pub iter_idxs: Vec<usize>,
    pub median_redshifts: Vec<f64>,
    pub distances: Vec<f64>,
    pub r50s: Vec<f64>,
    pub r100s: Vec<f64>,
    pub rsigmas: Vec<f64>,
    pub multiplicity: Vec<u32>,
    pub velocity_dispersion_gap: Vec<f64>,
    pub velocity_dispersion_gap_err: Vec<f64>,
    pub raw_masses: Vec<f64>,
    pub estimated_masses: Vec<f64>,
    pub bcg_idxs: Vec<usize>,
    pub bcg_ras: Vec<f64>,
    pub bcg_decs: Vec<f64>,
    pub bcg_redshifts: Vec<f64>,
    pub col_ras: Vec<f64>,
    pub col_decs: Vec<f64>,
    pub total_flux_proxies: Vec<f64>,
    pub total_absolute_mags: Vec<f64>,
}

/// Struct representing the pair catalog
pub struct PairCatalog {
    pub ids: Vec<i32>,
    pub idx_1: Vec<i32>,
    pub idx_2: Vec<i32>,
    pub projected_separation: Vec<f64>,
    pub velocity_separation: Vec<f64>,
    pub ra_bar: Vec<f64>,
    pub dec_bar: Vec<f64>,
    pub redshift_bar: Vec<f64>,
    pub total_absolute_mags: Vec<f64>,
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_gapper_velocity_basic() {
        let group = Group {
            redshift_members: vec![0.3, 0.2, 0.32, 0.5],
            velocity_errors: vec![50., 40., 30., 10.],
            ra_members: vec![0.2, 0.2, 0.2, 0.2],
            dec_members: vec![50., 50., 50., 50.],
            absolute_magnitude_members: vec![-18., -18., -18., -18.],
        };

        let result = group.velocity_dispersion_gapper();
        assert_eq!(result.0, 35908.750028805);
        assert_eq!(result.1, 5.70087712549569);
    }

    #[test]
    fn gapper_largeerror_0() {
        let group = Group {
            redshift_members: vec![0.3, 0.3, 0.3, 0.3],
            velocity_errors: vec![1000., 1000., 1000., 1000.],
            ra_members: vec![0.2, 0.2, 0.2, 0.2],
            dec_members: vec![50., 50., 50., 50.],
            absolute_magnitude_members: vec![-18., -18., -18., -18.],
        };
        let result = group.velocity_dispersion_gapper();
        assert_eq!(result.0, 0.0);
    }

    #[test]
    fn test_bright_central_galaxy_dominates_center() {
        let group = Group {
            ra_members: vec![180.0, 179.0, 181.0, 180.0],
            dec_members: vec![0.0, 1.0, -1.0, 0.0],
            redshift_members: vec![0.1; 4],
            absolute_magnitude_members: vec![-50.0, -18.0, -18.0, -18.0], // central very bright.
            velocity_errors: vec![50., 50., 50., 50.],
        };

        let iterative_center = group.calculate_iterative_center_idx();
        let (ra, dec) = (
            group.ra_members[iterative_center],
            group.dec_members[iterative_center],
        );
        // RA should be close to 180.0 and Dec to 0.0 (within ~0.01 deg)
        assert!((ra - 180.0).abs() < 1e-6, "RA deviated too far: {}", ra);
        assert!(dec.abs() < 1e-6, "Dec deviated too far: {}", dec);
    }

    #[test]
    fn testing_bcg_properties() {
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        // making a catalog with two groups, two binary pairs, and two single values.
        let catalog = GroupedGalaxyCatalog {
            ra: vec![
                23.1, 23., 23., 23., 43.1, 43., 43., 56.1, 56., 90., 90., 5., 270.,
            ],
            dec: vec![
                -20., -20.1, -20., -20., 0., 0.1, 0., 90., 90.1, 18., 18., -90., 56.,
            ],
            redshift: vec![
                0.21, 0.2, 0.2, 0.2, 0.41, 0.4, 0.4, 0.81, 0.8, 1.01, 1.0, 0.5, 0.5,
            ],
            absolute_magnitudes: vec![
                -23., -18., -18., -18., -23., -18., -18., -23., -18., -23., -18., -18., -18.,
            ],
            velocity_errors: vec![50.; 13],
            group_ids: vec![1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, -1, -1],
        };

        let results = catalog.calculate_group_properties(&cosmo);
        let ans = [0, 4, 7, 9];
        for (res, ans) in zip(results.bcg_idxs, ans) {
            assert_eq!(res, ans)
        }

        let result_ra = results.bcg_ras;
        let ans_ra = [23.1, 43.1, 56.1, 90.];

        let result_dec = results.bcg_decs;
        let ans_dec = [-20., 0., 90., 18.];

        let result_z = results.bcg_redshifts;
        let ans_redshift = [0.21, 0.41, 0.81, 1.01];

        for (res, ans) in zip(result_ra, ans_ra) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result_dec, ans_dec) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result_z, ans_redshift) {
            assert_eq!(res, ans)
        }
    }

    #[test]
    fn weird_group_different_to_master() {
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };
        let apparent_mags = [19.50136, 18.83319, 18.83163];
        let redshifts = vec![0.29022, 0.28753, 0.28734];

        let absolute_mags: Vec<f64> = apparent_mags
            .iter()
            .zip(redshifts.clone())
            .map(|(&mag, z)| mag - cosmo.distance_modulus(z))
            .collect();
        let group = Group {
            ra_members: vec![136.3741, 136.3661, 136.3610],
            dec_members: vec![1.698054, 1.719675, 1.711949],
            redshift_members: redshifts,
            absolute_magnitude_members: absolute_mags,
            velocity_errors: vec![50., 50., 50.],
        };

        let iterative_center = group.calculate_iterative_center_idx();
        let (ra, dec) = (
            group.ra_members[iterative_center],
            group.dec_members[iterative_center],
        );
        assert_eq!(ra, group.ra_members[1]);
        assert_eq!(dec, group.dec_members[1]);
    }

    #[test]
    fn testing_flux_vals() {
        let group = Group {
            ra_members: vec![23., 23.2, 22.9, 24.0],
            dec_members: vec![-23., -23.2, -23.2, -23.],
            redshift_members: vec![0.2, 0.2, 0.2, 0.2],
            absolute_magnitude_members: vec![-18., -18., -18., -18.],
            velocity_errors: vec![50., 50., 50., 50.],
        };

        let result_proxy = group.flux_proxy();
        let result_mag = group.total_absolute_magnitude();
        let answer_proxy = 4677997564.079487;
        let answer_mag = -19.505149978319906;
        assert_eq!(result_proxy, answer_proxy);
        assert_eq!(result_mag, answer_mag);
    }

    #[test]
    fn comparing_radii_to_r() {
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        let group = Group {
            ra_members: vec![23., 23.2, 22.9, 24.0],
            dec_members: vec![-23., -23.2, -23.2, -23.],
            redshift_members: vec![0.2, 0.2, 0.2, 0.2],
            absolute_magnitude_members: vec![-18., -18., -18., -18.],
            velocity_errors: vec![50., 50., 50., 50.],
        };

        let group_ra = 23.;
        let group_dec = -23.;
        let group_z = 0.2;

        let answer = [3.505739, 4.243428, 13.121179];
        let result = group.calculate_radius(group_ra, group_dec, group_z, &cosmo);
        for (r, a) in zip(result, answer) {
            assert!((r - a).abs() < 1e-6)
        }
    }

    #[test]
    fn test_group_catalog() {
        let cosmo = Cosmology {
            omega_m: 0.3,
            omega_k: 0.,
            omega_l: 0.7,
            h0: 70.,
        };

        // making a catalog with two groups, a binary pair, and two single values.
        let catalog = GroupedGalaxyCatalog {
            ra: vec![
                23., 23., 23., 23., 43., 43., 43., 56., 56., 90., 90., 5., 270.,
            ],
            dec: vec![
                -20., -20., -20., -20., 0., 0., 0., 90., 90., 18., 18., -90., 56.,
            ],
            redshift: vec![
                0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.8, 0.8, 1.0, 1.0, 0.5, 0.5,
            ],
            absolute_magnitudes: vec![-18.; 13],
            velocity_errors: vec![50.; 13],
            group_ids: vec![1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, -1, -1],
        };

        let result = catalog.calculate_group_properties(&cosmo);
        for (res_dec, ans_dec) in zip(result.iter_decs, vec![-20., 0., 90., 18.]) {
            assert!((res_dec - ans_dec).abs() < 1e-8);
        }

        for (res_ra, ans_ra) in zip(result.iter_ras, vec![23., 43., 56., 90.]) {
            assert!((res_ra - ans_ra).abs() < 1e-8);
        }

        for (res_red, ans_red) in zip(&result.median_redshifts, vec![0.2, 0.4, 0.8, 1.0]) {
            assert!((res_red - ans_red).abs() < 1e-8);
        }

        for (res_mult, ans_mult) in zip(result.multiplicity, vec![4, 3, 2, 2]) {
            assert_eq!(res_mult, ans_mult);
        }

        for (res_ids, ans_ids) in zip(result.ids, vec![1, 2, 3]) {
            assert_eq!(res_ids, ans_ids)
        }
        let ans_distances: Vec<f64> = result
            .median_redshifts
            .into_iter()
            .map(|z| cosmo.comoving_distance(z))
            .collect();

        for (res_distance, ans_distance) in zip(result.distances, ans_distances) {
            assert_eq!(res_distance, ans_distance)
        }

        for (res_50, ans_50) in zip(result.r50s, vec![0., 0., 0.]) {
            assert!((res_50 - ans_50).abs() < 1e-7)
        }

        for (res_100, ans_100) in zip(result.r100s, vec![0., 0., 0.]) {
            assert!((res_100 - ans_100).abs() < 1e-7)
        }

        for (res_sig, ans_sig) in zip(result.rsigmas, vec![0., 0., 0.]) {
            assert!((res_sig - ans_sig).abs() < 1e-7)
        }
    }

    #[test]
    fn testing_pair_catalog() {
        // making a catalog with two groups, two binary pairs, and two single values.
        let catalog = GroupedGalaxyCatalog {
            ra: vec![
                23., 23., 23., 23., 43., 43., 43., 56., 56., 90., 90., 5., 270.,
            ],
            dec: vec![
                -20., -20., -20., -20., 0., 0., 0., 90., 90., 18., 18., -90., 56.,
            ],
            redshift: vec![
                0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.8, 0.81, 1.1, 1.0, 0.5, 0.5,
            ],
            absolute_magnitudes: vec![-18.; 13],
            velocity_errors: vec![50.; 13],
            group_ids: vec![1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, -1, -1],
        };
        let result = catalog.calculate_pair_properties();
        let answer_id = [3, 4];
        let answer_id_1 = [7, 9];
        let answer_id_2 = [8, 10];
        let answer_projected_separation = [0., 0.];
        let answer_vel_sep = [0.01, 0.1];

        for (res, ans) in zip(result.ids, answer_id) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result.idx_1, answer_id_1) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result.idx_2, answer_id_2) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result.projected_separation, answer_projected_separation) {
            assert_eq!(res, ans)
        }
        for (res, ans) in zip(result.velocity_separation, answer_vel_sep) {
            assert!((res - ans).abs() < 1e-7)
        }
    }
}
