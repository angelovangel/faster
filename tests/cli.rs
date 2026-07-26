use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs


#[test]
fn find_content_in_table() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-ts").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("10\t18931\t0\t165"));

    Ok(())
}

#[test]
fn find_content_in_len() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-l").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("249"));

    Ok(())
}

// --- filterl ---

#[test]
fn filterl_keep_longer_than_1000() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-f").arg("1000").arg("tests/test.fastq");

    // 3 reads are longer than 1000 bp (2440, 5263, 8490)
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("3 kept"));

    Ok(())
}

#[test]
fn filterl_keep_shorter_than_500() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-f").arg("-500").arg("tests/test.fastq");

    // 6 reads are shorter than 500 bp (165, 240, 249, 391, 421, 453)
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("6 kept"));

    Ok(())
}

// --- filterq ---

#[test]
fn filterq_keep_above_q10() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-w").arg("10").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stderr(predicate::str::contains("6 kept"));

    Ok(())
}

// --- dust ---

#[test]
fn dust_filter_default_threshold() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-d").arg("2.0").arg("tests/test.fastq");

    // real biological reads should all pass the dust filter
    cmd.assert()
        .success()
        .stderr(predicate::str::contains("10 kept"));

    Ok(())
}

// --- nx ---

#[test]
fn nx_n50() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-x").arg("0.5").arg("tests/test.fastq");

    // N50 of the test file is 5263
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("N50\t5263"));

    Ok(())
}

// --- qyield ---

#[test]
fn qyield_q20() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-y").arg("20").arg("tests/test.fastq");

    // output should start with "Q20" and contain a percentage
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Q20"));

    Ok(())
}
