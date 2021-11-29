use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::process::Command; // Run programs


#[test]

fn find_content_in_output() -> Result<(), Box<dyn std::error::Error>> {

    let mut cmd = Command::cargo_bin("faster")?;
    cmd.arg("-ts").arg("tests/test.fastq");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains( "11313\t60245316\t0\t143") );

    Ok(())
}