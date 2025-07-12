use fof::group_properties::GroupedGalaxyCatalog;
use fof::cosmology_funcs::Cosmology;



#[test]
fn test_group_properties() {
    let rdr = csv::Reader::from_path("tests/test_group_properties.csv");
    
    let mut ras = Vec::new();
    let mut decs = Vec::new();
    let mut redshifts = Vec::new();
    let mut ab_mags = Vec::new();
    let mut vel_errs = Vec::new();
    let mut group_ids = Vec::new();

    for result in rdr.expect("No file").records() {
        let record = result.unwrap();
        ras.push(record[1].parse::<f64>().unwrap());
        decs.push(record[2].parse::<f64>().unwrap());
        redshifts.push(record[3].parse::<f64>().unwrap());
        ab_mags.push(record[4].parse::<f64>().unwrap());
        vel_errs.push(record[5].parse::<f64>().unwrap());
        group_ids.push(record[6].parse::<i32>().unwrap());
    }

    let catalog = GroupedGalaxyCatalog {
        ra: ras,
        dec: decs,
        redshift: redshifts,
        absolute_magnitudes: ab_mags,
        velocity_errors: vel_errs,
        group_ids,
    };

    let cosmo = &Cosmology {
        omega_m: 0.3,
        omega_k: 0.0,
        omega_l: 0.7,
        h0: 0.7,
    };
    let group_catalog = catalog.calculate_group_properties(cosmo);

    assert!(group_catalog.multiplicity[0] > 0, "No groups were found");

} 