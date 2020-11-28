//use std::io;
use std::io::BufReader;
use std::fs;
use std::process;
use flate2::bufread;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use bio::seq_analysis::gc::gc_content;

extern crate clap;
use clap::{Arg, App};

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

// just average
fn mean(numbers: &[i32]) -> f32 {
    numbers.iter().sum::<i32>() as f32 / numbers.len() as f32
}

// median is not precise for even numbers, but ok for this application, keep it as i32
fn median(numbers: &mut [i32]) -> i32 {
    numbers.sort();
    let mid = numbers.len() / 2;
    numbers[mid]
}

// n50 , TODO - make it nX with a second parameter
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

// get number of bases with q >= value
fn get_qual_bases(q: &[u8], qx: u8) -> i32 {
    let mut n = 0;
    for &item in q.iter()
     {
        if *&item >= qx {
            n += 1
        }
    }
    n
}

fn main() {

    let matches = App::new("fastq-stats")
                        .version("0.1.3")
                        .author("Angel Angelov <aangeloo@gmail.com>")
                        .about("fast statistics and more for 1 fastq file")
                        .arg(Arg::with_name("stats")
                            .short("s")
                            .long("stats")
                            .conflicts_with_all(&["len", "gc"]) //conflicts_with is two-way, so no need for the others
                            .help("Output a table with statistics about the fastq file"))
                        .arg(Arg::with_name("len")
                            .short("l")
                            .long("len")
                            .conflicts_with("gc")
                            .help("Output read lengths, one line per read"))
                        .arg(Arg::with_name("gc")
                            .short("g")
                            .long("gc")
                            .help("Output gc values, one line per read"))
                        .arg(Arg::with_name("INPUT")
                            .help("path to fastq file")
                            .required(true)
                            .index(1))
                        .get_matches();
    //println!("Working on {}", matches.value_of("INPUT").unwrap());
    // read file
    // Here we can call .unwrap() because the argument is required.
    let infile = matches.value_of("INPUT").unwrap().to_string();

    let mut reader = fastq::Reader::new(get_fastq_reader(&infile));
    let mut record = fastq::Record::new();

    let mut reads = 0;
    let mut bases = 0;
    let mut qual20 = 0;
    let mut qual30 = 0;
    let mut minlen: i32 = i32::MAX;
    let mut maxlen = 0;
    let mut len_vector: Vec<i32> = Vec::new();

    // here discriminate output based on arguments
    //case len

    if matches.is_present("len") {
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
        process::exit(0);

    // case gc
    } else if matches.is_present("gc") {
        reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");
        while !record.is_empty() {
            let seq = record.seq();
            println!("{}", gc_content(seq));

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    }

    // normal case, output table
    reader
        .read(&mut record)
        .expect("Failed to parse fastq record!");
    while !record.is_empty() {
        //let seq = record.seq();
        let len = record.seq().len() as i32;
        
        reads += 1;
        bases += len;
        qual20 += get_qual_bases(record.qual(), 53); // 33 offset
        qual30 += get_qual_bases(record.qual(), 63);
        minlen = len.min(minlen);
        maxlen = len.max(maxlen);
        len_vector.push(len);

        reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");
    }

    let mean_len = mean(&len_vector);
    let median_len = median(&mut len_vector);
    let n50 = n50(&mut len_vector, 0.5); // use 0.1 for N90!!!
    let q20 = qual20 as f64 / bases as f64 * 100.0;
    let q30 = qual30 as f64 / bases as f64 * 100.0;

    println!("file\treads\tbases\tmin_len\tmax_len\tmean_len\tmedian_len\tN50\tQ20_percent\tQ30_percent");
    println!("{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{:.2}\t{:.2}", infile, reads, bases, minlen, maxlen, mean_len, median_len, n50, q20, q30);
}
