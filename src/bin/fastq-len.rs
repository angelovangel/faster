// get len per read in fastq file
use std::io::BufReader;
use std::env;
use std::fs;
use std::process;
use flate2::bufread;
use bio::io::fastq;
use bio::io::fastq::FastqRead;

// fastq reader, file as arg, decide based on extension
fn get_fastq_reader(path: &String) -> Box<dyn (::std::io::Read)> {
    if path.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
    } else {
        let f = fs::File::open(path).unwrap();
        Box::new(BufReader::new(f))
        //Box::new(fs::File::open(path).unwrap())
    }
}


fn main() {

    // read file
    let args: Vec<String> = env::args().collect();
    let arg1 = &args[1];
    
    if arg1.starts_with("--help") {
        println!("Usage: fastq-gc file.fastq[.gz] \n \n\
                  Input file may be compressed using gzip, \n\
                  in which case it must end with '.gz'");
        process::exit(1);
    }

    let mut reader = fastq::Reader::new(get_fastq_reader(arg1));
    let mut record = fastq::Record::new();

    reader
        .read(&mut record)
        .expect("Failed to parse fastq record!");
    
        while !record.is_empty() {
        let len = record.seq().len() as i32;
        println!("{}", len);

        reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");
    }
    
    
}