use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs


#[test]

fn find_content_in_table() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-ts").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains( "10\t18931\t0\t165") );

    Ok(())
}

#[test]
fn find_content_in_len() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-l").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains( "249") );

    Ok(())
}
