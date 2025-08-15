use std::{time::Instant};

use fof::link_finder::{find_links};

#[test]
fn recovering_ffl1() {
    // read in the data
    let rdr = csv::Reader::from_path("tests/test_links.csv");

    let mut ras = Vec::new();
    let mut decs = Vec::new();
    let mut distances = Vec::new();
    let mut pos = Vec::new();
    let mut los: Vec<f64> = Vec::new();

    for result in rdr.expect("No file").records() {
        let record = result.unwrap();
        ras.push(record[0].parse::<f64>().unwrap());
        decs.push(record[1].parse::<f64>().unwrap());
        distances.push(record[2].parse::<f64>().unwrap());
        pos.push(record[3].parse::<f64>().unwrap());
        los.push(record[4].parse::<f64>().unwrap());
    }

    let now_new = Instant::now();
    let _answer = find_links(
        ras.clone(),
        decs.clone(),
        distances.clone(),
        pos.clone(),
        los.clone(),
    );
    let tree_time = now_new.elapsed();


    println!("Time: {:?}", tree_time);
}
