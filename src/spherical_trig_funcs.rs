use libm::{asin, atan2};

// assuming a unit sphere
pub fn convert_equitorial_to_cartesian(ra_deg: &f64, dec_deg: &f64) -> [f64; 3] {
    let ra_radians = ra_deg.to_radians();
    let dec_radians = dec_deg.to_radians();
    let x = dec_radians.cos() * ra_radians.cos();
    let y = dec_radians.cos() * ra_radians.sin();
    let z = dec_radians.sin();
    [x, y, z]
}

// assuming a unit sphere
pub fn convert_cartesian_to_equitorial(x: &f64, y: &f64, z: &f64) -> [f64; 2] {
    let radius = (x.powi(2) + y.powi(2) + z.powi(2)).sqrt();
    let ra_radian = atan2(*y, *x);
    let dec_radian = asin(z / radius);
    [ra_radian.to_degrees(), dec_radian.to_degrees()]
}

// Converts RA, Dec (degrees) and comoving distance to 3D Cartesian coordinates.
pub fn convert_equitorial_to_cartesian_scaled(
    ra_deg: f64,
    dec_deg: f64,
    distance: f64,
) -> [f64; 3] {
    let ra_rad = ra_deg.to_radians();
    let dec_rad = dec_deg.to_radians();
    let x = distance * dec_rad.cos() * ra_rad.cos();
    let y = distance * dec_rad.cos() * ra_rad.sin();
    let z = distance * dec_rad.sin();
    [x, y, z]
}

// The three dimensional euclidean distance
pub fn euclidean_distance_3d(point_1: &[f64; 3], point_2: &[f64; 3]) -> f64 {
    ((point_1[0] - point_2[0]).powi(2)
        + (point_1[1] - point_2[1]).powi(2)
        + (point_1[2] - point_2[2]).powi(2))
    .sqrt()
}

/// Chord length. Given the angular separation.
pub fn chord_distance(angular_separation_degrees: f64) -> f64 {
    assert!((0. ..=180.).contains(&angular_separation_degrees));
    2. * (angular_separation_degrees.to_radians() / 2.).sin()
}

/// projected separation between two galaxies. Angular scale only
pub fn calculate_projected_separation(ra_deg_1: &f64, dec_deg_1: &f64, ra_deg_2: &f64, dec_deg_2: &f64) -> f64 {
    let point_1 = convert_equitorial_to_cartesian(ra_deg_1, dec_deg_1);
    let point_2 = convert_equitorial_to_cartesian(ra_deg_2, dec_deg_2);
    euclidean_distance_3d(&point_1, &point_2)
}

#[cfg(test)]
mod test {
    use super::*;
    use std::iter::zip;

    const EPSILON: f64 = 1e-6;

    #[test]
    fn test_equatorial_to_cartesian_unit_sphere() {
        let ra = 0.0;
        let dec = 0.0;
        let result = convert_equitorial_to_cartesian(&ra, &dec);
        assert!((result[0] - 1.0).abs() < EPSILON);
        assert!((result[1]).abs() < EPSILON);
        assert!((result[2]).abs() < EPSILON);
    }

    #[test]
    fn test_cartesian_to_equatorial_unit_sphere() {
        let x = 1.0;
        let y = 0.0;
        let z = 0.0;
        let result = convert_cartesian_to_equitorial(&x, &y, &z);
        assert!((result[0] - 0.0).abs() < EPSILON); // RA
        assert!((result[1] - 0.0).abs() < EPSILON); // Dec
    }

    #[test]
    fn test_round_trip_conversion() {
        let ra = 123.4;
        let dec = -56.78;
        let cart = convert_equitorial_to_cartesian(&ra, &dec);
        let equi = convert_cartesian_to_equitorial(&cart[0], &cart[1], &cart[2]);

        // RA can wrap around so we normalize to [0, 360)
        let mut ra_norm = equi[0];
        if ra_norm < 0.0 {
            ra_norm += 360.0;
        }

        assert!((ra_norm - ra).abs() < EPSILON);
        assert!((equi[1] - dec).abs() < EPSILON);
    }

    #[test]
    fn test_equatorial_to_cartesian_scaled() {
        let ra = [120., 122., 124.];
        let dec = [-56., -34., -30.];
        let distance = [10., 11., 12.];
        let result_1 = convert_equitorial_to_cartesian_scaled(ra[0], dec[0], distance[0]);
        let result_2 = convert_equitorial_to_cartesian_scaled(ra[1], dec[1], distance[1]);
        let result_3 = convert_equitorial_to_cartesian_scaled(ra[2], dec[2], distance[2]);

        let answer_1 = [-2.795965, 4.842753, -8.290376];
        let answer_2 = [-4.832553, 7.733701, -6.151122];
        let answer_3 = [-5.811303, 8.615611, -6.000000];

        for (r, a) in zip(result_1, answer_1) {
            assert!((r - a).abs() < EPSILON);
        }

        for (r, a) in zip(result_2, answer_2) {
            assert!((r - a).abs() < EPSILON);
        }

        for (r, a) in zip(result_3, answer_3) {
            assert!((r - a).abs() < EPSILON);
        }
    }

    #[test]
    fn test_euclidean_distance_3d() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 0.0, 0.0];
        let dist = euclidean_distance_3d(&a, &b);
        assert!((dist - 1.0).abs() < EPSILON);

        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 6.0, 3.0];
        let dist = euclidean_distance_3d(&a, &b);
        assert!((dist - 5.0).abs() < EPSILON);
    }

    #[test]
    fn test_chord_distance() {
        let angular_separation = [0., 20., 30., 180.];
        let results = angular_separation.iter().map(|&a| chord_distance(a));
        let answers = [0., 0.347296355, 0.517638090, 2.];
        for (res, ans) in zip(results, answers) {
            assert!((res - ans).abs() < 1e-7)
        }
    }
}
