use std::{iter::zip, time::Instant};

use fof::link_finder::{ffl1, find_links};

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
    let answer = find_links(
        ras.clone(),
        decs.clone(),
        distances.clone(),
        pos.clone(),
        los.clone(),
    );
    let tree_time = now_new.elapsed();

    let now_old = Instant::now();
    let result = ffl1(
        ras.clone(),
        decs.clone(),
        distances.clone(),
        pos.clone(),
        los.clone(),
    );
    let old_time = now_old.elapsed();

    let now_class = Instant::now();
    //let classic = ffl1(ras.clone(), decs.clone(), distances.clone(), pos.clone(), los.clone());
    let classic_time = now_class.elapsed();

    println!("Old time: {:?}", old_time);
    println!("new time: {:?}", tree_time);
    println!("classic time: {:?}", classic_time);

    let mut answer_left: Vec<usize> = answer.iter().map(|u| u.0).collect();
    let mut result_left: Vec<usize> = result.iter().map(|u| u.0).collect();
    let mut answer_right: Vec<usize> = answer.iter().map(|u| u.1).collect();
    let mut result_right: Vec<usize> = result.iter().map(|u| u.1).collect();

    answer_left.sort();
    result_left.sort();
    answer_right.sort();
    result_right.sort();

    for (res, ans) in zip(result_left.clone(), answer_left.clone()) {
        assert_eq!(res, ans)
    }

    for (res, ans) in zip(result_right.clone(), answer_right.clone()) {
        assert_eq!(res, ans)
    }
}
