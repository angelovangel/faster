use std::{io, io::BufReader, fs, process};
use flate2::bufread;
use bio::io::{fastq, fastq::FastqRead};
use bio::seq_analysis::gc::gc_content;
use rand::seq::IteratorRandom;

extern crate clap;
use clap::{Arg, App, ArgGroup};
// own functions
mod modules;

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

// subsample records (reads) from a fastq file, given number of reads
fn samplefq(path: &String, n: usize) {

    let records = fastq::Reader::new(get_fastq_reader(path)).records();
    let mut writer = fastq::Writer::new(io::stdout());
    let mut rng = rand::thread_rng();

    // sampled is a vector of result types
    let sampled = records.choose_multiple(&mut rng, n);
    for record in sampled {
        // get the record
        let r = record.unwrap();
        writer
            .write_record(&r)
            .expect("Error in writing fastq record!");
    }
}

fn main() {

    let matches = App::new("faster")
                        .version("0.1.3")
                        .author("Angel Angelov <aangeloo@gmail.com>")
                        .about("fast statistics and more for 1 fastq file")

                        .arg(Arg::with_name("table")
                            .short("t")
                            .long("table")
                            .help("Output a (tab separated) table with statistics"))

                        .arg(Arg::with_name("len")
                            .short("l")
                            .long("len")
                            .help("Output read lengths, one line per read"))

                        .arg(Arg::with_name("gc")
                            .short("g")
                            .long("gc")
                            .help("Output gc values, one line per read"))

                        .arg(Arg::with_name("qscore")
                            .short("q")
                            .long("qscore")
                            .help("Output 'mean' read qscores, one line per read. For this, the mean of the base probabilities for each read is calculated, and the result is converted back to a phred score"))

                        .arg(Arg::with_name("sample")
                        .short("s")
                        .long("sample")
                        .takes_value(true)
                        .help("Sub-sample sequences by proportion"))

                        .arg(Arg::with_name("filter")
                            .short("f")
                            .long("filter")
                            .takes_value(true)
                            .help("Filter reads based on length - only reads with length greater than [integer] are written to stdout"))

                        .arg(Arg::with_name("INPUT")
                            .help("path to fastq file")
                            .required(true)
                            .index(1))

                        // this group makes one and only one arg from the set required, avoid defining conflicts_with
                        .group(ArgGroup::with_name("group")
                        .required(true).args(&["table", "len", "gc", "qscore", "filter", "sample"]))
                        .get_matches();
    //println!("Working on {}", matches.value_of("INPUT").unwrap());
    // read file
    // Here we can call .unwrap() because the argument is required.
    let infile = matches.value_of("INPUT").unwrap().to_string();

    let mut reader = fastq::Reader::new(get_fastq_reader(&infile));
    let mut record = fastq::Record::new();
    let mut writer = fastq::Writer::new(io::stdout());

    // read first record, then check for arguments
    reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");

    // stay in loop until all records are read
    //case len
    if matches.is_present("len") {

        while !record.is_empty() {
            let len = record.seq().len() as i64;
            println!("{}", len);

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);

    // case gc
    } else if matches.is_present("gc") {

        while !record.is_empty() {
            let seq = record.seq();
            println!("{}", gc_content(seq));

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);

    // case qscore
    } else if matches.is_present("qscore") {

        while !record.is_empty() {
            let qscore = modules::qscore_probs(record.qual()) / record.seq().len() as f32;
            println!("{:.4}", -10.0 * qscore.log10());

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    
    // case filter
    } else if matches.is_present("filter") {
        let filterlen = matches.value_of("filter").unwrap();
        // why so much code to parse an integer from arg!
        let input_opt = filterlen.trim().parse::<i32>();
        let input_int = match input_opt {
            Ok(input_int) => input_int,
            Err(e) => {
                println!("Please use an integer for -f ({})", e);
                return;
            }
        };

        while !record.is_empty() {
            let seqlen = record.seq().len() as i32;
            // try branchless?
            if seqlen > input_int {
                writer
                    .write_record(&record)
                    .expect("Failed to write file!");
            }

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
        
    } else if matches.is_present("sample") {
        // get n reads first
        let mut reads: i64 = 0;
        while !record.is_empty() {
            reads += 1;
            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }

        // parse fraction
        let fraction = matches
            .value_of("sample")
            .unwrap()
            .trim()
            .parse::<f32>()
            .unwrap();
        // how many reads to sample? also check if fraction is 0..1
        match fraction {
            // see <https://stackoverflow.com/a/58434531/8040734>
            x if (0.0..1.0).contains(&x) => {
                let nreads = reads as f32 * fraction;
                samplefq(&infile, nreads as usize);
                //println!("nreads: {}", nreads as usize);
                process::exit(0);
            },
            _ => eprintln!("The subsample fraction should be between 0.0 and 1.0!")
        }
        
    }
        
    // normal case, output table
    let mut reads: i64 = 0;
    let mut bases: i64 = 0;
    let mut qual20: i64 = 0;
    let mut qual30: i64 = 0;
    let mut minlen: i64 = i64::MAX;
    let mut maxlen = 0;
    let mut len_vector: Vec<i64> = Vec::new();

    while !record.is_empty() {
        //let seq = record.seq();
        let len = record.seq().len() as i64; // here have to accomodate bigger numbers, as bases can get > 2^32
        
        reads += 1;
        bases += len;
        qual20 += modules::get_qual_bases(record.qual(), 53); // 33 offset
        qual30 += modules::get_qual_bases(record.qual(), 63);
        minlen = len.min(minlen);
        maxlen = len.max(maxlen);
        len_vector.push(len);

        reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");
    }

    let mean_len = modules::mean(&len_vector);
    let quart1 = modules::quartiles(&mut len_vector, 1);
    let quart2 = modules::quartiles(&mut len_vector, 2);
    let quart3 = modules::quartiles(&mut len_vector, 3);
    let n50 = modules::get_nx(&mut len_vector, 0.5); // use 0.1 for N90!!!
    let q20 = qual20 as f64 / bases as f64 * 100.0;
    let q30 = qual30 as f64 / bases as f64 * 100.0;

    println!("file\treads\tbases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent");
    println!("{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}", infile, reads, bases, minlen, maxlen, mean_len, quart1, quart2, quart3, n50, q20, q30);
}
// END