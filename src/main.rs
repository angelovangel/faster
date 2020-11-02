//use std::io;
use std::io::BufReader;
use std::env;
use std::fs;
use std::process;
use flate2::bufread;
use bio::io::fastq;
use bio::io::fastq::FastqRead;

// fastq reader, file as arg, decide based on extension
fn get_fastq_reader(path: &String) -> Box<::std::io::Read> {
    if path.ends_with(".gz") {
        let f = fs::File::open(path).unwrap();
        Box::new(bufread::MultiGzDecoder::new(BufReader::new(f)))
    } else {
        Box::new(fs::File::open(path).unwrap())
    }
}
// take a fastq on stdin and output a tab separated table with useful info

// just average
fn average(numbers: &[i32]) -> f32 {
    numbers.iter().sum::<i32>() as f32 / numbers.len() as f32
}

// median is not precise for even numbers, but ok for this application, keep it as i32
fn median(numbers: &mut [i32]) -> i32 {
    numbers.sort();
    let mid = numbers.len() / 2;
    numbers[mid]
}

// n5o , TODO - make it nX with a second parameter
fn n50(numbers: &mut [i32], fraction: f32) -> i32 {

    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.iter().sum::<i32>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers.iter()
        .scan(0, |sum, i | { *sum += i; Some(*sum) })
        .collect::<Vec<_>>();
    let n50_index = cumsum
        .iter()
        .position(|&x| x > halfsum as i32)
        .unwrap();

    numbers[n50_index]
}

fn main() {

    // read file
    let args: Vec<String> = env::args().collect();
    let arg1 = &args[1];

    if arg1.starts_with("--help") {
        println!("Usage: fastq-stats fastq-file.fastq[.gz] \n \n\
                  Input files may be compressed using gzip, \n\
                  in which case they must end with '.gz'");
        process::exit(0);
    }

    let mut reader = fastq::Reader::new(get_fastq_reader(arg1));
    let mut record = fastq::Record::new();

    let mut reads = 0;
    let mut bases = 0;
    let mut minlen: i32 = i32::MAX;
    let mut maxlen = 0;
    let mut len_vector: Vec<i32> = Vec::new();


    reader
        .read(&mut record)
        .expect("Failed to parse fastq record!");

    while !record.is_empty() {
        let seq = record.seq();
        let len = record.seq().len() as i32;

        reads += 1;
        bases += seq.len();
        minlen = len.min(minlen);
        maxlen = len.max(maxlen);
        len_vector.push(len);

        reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");
    }

    let av_len = average(&len_vector);
    let median_len = median(&mut len_vector);
    let n50 = n50(&mut len_vector, 0.5); // use 0.1 for N90!!!

    println!("reads\tbases\tminlen\tav_len\tmedian_len\tmaxlen\tN50");
    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}", reads, bases, minlen, av_len, median_len, maxlen, n50);
}

