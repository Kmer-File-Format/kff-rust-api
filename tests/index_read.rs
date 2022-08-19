//! Check example produce the good output

use std::io::Read;
use std::process::{Command, Stdio};

#[test]
#[cfg(not(tarpaulin))]
fn get_index() -> kff::error::Result<()> {
    let args = vec![
        "run",
        "--example",
        "index_read",
        "--",
        "-i",
        "tests/data/index.kff",
        "-s",
    ];

    let mut child = Command::new("cargo")
        .args(&args)
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create cargo example child");

    if !child.wait().expect("Error durring yacrd run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout)?;
        child.stderr.unwrap().read_to_string(&mut stderr)?;

        panic!("\nstdout: {}\nstderr: {}", stdout, stderr);
    }

    let mut value_str = String::new();
    child.stdout.unwrap().read_to_string(&mut value_str)?;
    value_str.pop().unwrap();

    let value: Vec<Vec<u8>> = value_str
        .split('\n')
        .map(|x| x.bytes().collect::<Vec<u8>>())
        .collect();

    assert_eq!(
        value,
        vec![
            b"index, type, offset".to_vec(),
            b"0,v,23".to_vec(),
            b"1,m,82".to_vec(),
            b"2,m,97".to_vec(),
            b"3,m,113".to_vec(),
            b"4,m,130".to_vec(),
            b"5,m,145".to_vec(),
            b"6,m,161".to_vec(),
            b"7,m,177".to_vec(),
            b"8,m,192".to_vec(),
            b"9,m,208".to_vec(),
            b"10,m,223".to_vec(),
            b"11,m,240".to_vec(),
            b"12,m,256".to_vec(),
            b"13,m,285".to_vec(),
            b"14,m,300".to_vec(),
            b"15,m,316".to_vec(),
            b"16,m,333".to_vec(),
            b"17,m,348".to_vec(),
            b"18,m,363".to_vec(),
            b"19,m,380".to_vec(),
            b"20,m,396".to_vec(),
            b"21,m,412".to_vec(),
            b"22,m,427".to_vec(),
            b"23,m,442".to_vec(),
            b"24,m,459".to_vec(),
            b"25,m,476".to_vec(),
            b"26,m,493".to_vec(),
            b"27,m,508".to_vec(),
            b"28,m,525".to_vec(),
            b"29,m,541".to_vec(),
            b"30,m,556".to_vec(),
            b"31,m,572".to_vec(),
            b"32,m,601".to_vec(),
            b"33,m,618".to_vec(),
            b"34,m,635".to_vec(),
            b"35,m,650".to_vec(),
            b"36,m,665".to_vec(),
            b"37,m,680".to_vec(),
            b"38,m,706".to_vec(),
            b"39,m,732".to_vec(),
            b"40,m,749".to_vec(),
            b"41,m,765".to_vec(),
            b"42,m,782".to_vec(),
            b"43,m,798".to_vec(),
            b"44,m,815".to_vec(),
            b"45,m,830".to_vec(),
            b"46,m,852".to_vec(),
            b"47,m,869".to_vec(),
            b"48,m,890".to_vec(),
            b"49,m,905".to_vec(),
        ]
    );

    Ok(())
}

#[test]
#[cfg(not(tarpaulin))]
fn get_section() -> kff::error::Result<()> {
    let args = vec![
        "run",
        "--example",
        "index_read",
        "--",
        "-i",
        "tests/data/index.kff",
        "-t",
        "5",
    ];

    let mut child = Command::new("cargo")
        .args(&args)
        .stderr(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("Couldn't create cargo example child");

    if !child.wait().expect("Error durring yacrd run").success() {
        let mut stdout = String::new();
        let mut stderr = String::new();

        child.stdout.unwrap().read_to_string(&mut stdout)?;
        child.stderr.unwrap().read_to_string(&mut stderr)?;

        panic!("\nstdout: {}\nstderr: {}", stdout, stderr);
    }

    let mut value_str = String::new();
    child.stdout.unwrap().read_to_string(&mut value_str)?;
    value_str.pop().unwrap();

    let value: Vec<Vec<u8>> = value_str
        .split('\n')
        .map(|x| x.bytes().collect::<Vec<u8>>())
        .collect();

    assert_eq!(
        value,
        vec![
            b"GAGCCTTAACGTAGC".to_vec(),
            b"AGCCTTAACGTAGCG".to_vec(),
            b"GCCTTAACGTAGCGC".to_vec(),
            b"CCTTAACGTAGCGCC".to_vec(),
            b"CTTAACGTAGCGCCA".to_vec(),
        ]
    );

    Ok(())
}
