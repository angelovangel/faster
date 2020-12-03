use std::io;
use std::io::BufReader;
use std::fs;
use std::process;
use flate2::bufread;
use bio::io::fastq;
use bio::io::fastq::FastqRead;
use bio::seq_analysis::gc::gc_content;

extern crate clap;
use clap::{Arg, App, ArgGroup};

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
fn mean(numbers: &[i64]) -> f64 {
    numbers.iter().sum::<i64>() as f64 / numbers.len() as f64
}

// median is not precise for vectors with even numbers, but ok for this application, keep it as i32
fn quartiles(numbers: &mut [i64], q: i8) -> i64 {
    numbers.sort_unstable();
    match q {
        1 => {
            let index = numbers.len() / 4;
            return numbers[index]
        },
        2 => {
            let index = numbers.len() / 2;
            return numbers[index]
        },
        3 => {
            // avoid having to use f64
            let index1 = numbers.len() / 4;
            let index2 = numbers.len() / 2;
            return numbers[index1 + index2]
        },
        _ => 42 //:)
    }
    // first quartile
    
    
}

// n50 , TODO - make it nX with a second parameter
fn n50(numbers: &mut [i64], fraction: f32) -> i64 {

    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers.iter()
        .scan(0, |sum, i | { *sum += i; Some(*sum) })
        .collect::<Vec<_>>();
    let n50_index = cumsum
        .iter()
        .position(|&x| x > halfsum as i64)
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

// to get mean of q scores from a record - first convert to prob, calc mean, then back to phred
// this fn reads phred and converts to probs and returns their sum
// 
fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum= 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred/10.0);
        qprob_sum += prob
    }
    qprob_sum
}

fn main() {

    let matches = App::new("faster")
                        .version("0.1.3")
                        .author("Angel Angelov <aangeloo@gmail.com>")
                        .about("fast statistics and more for 1 fastq file")

                        .arg(Arg::with_name("stats")
                            .short("s")
                            .long("stats")
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
                        .required(true).args(&["stats", "len", "gc", "qscore", "filter"]))
                        .get_matches();
    //println!("Working on {}", matches.value_of("INPUT").unwrap());
    // read file
    // Here we can call .unwrap() because the argument is required.
    let infile = matches.value_of("INPUT").unwrap().to_string();

    let mut reader = fastq::Reader::new(get_fastq_reader(&infile));
    let mut record = fastq::Record::new();
    let mut writer = fastq::Writer::new(io::stdout());


    reader
            .read(&mut record)
            .expect("Failed to parse fastq record!");

    // here discriminate output based on arguments
    //case len
    if matches.is_present("len") {

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
            let qscore = qscore_probs(record.qual()) / record.seq().len() as f32;
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
                writer.write_record(&record)
                    .expect("Failed to write file!");
            }

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    }
        
    // normal case, output table
    let mut reads = 0;
    let mut bases = 0;
    let mut qual20 = 0;
    let mut qual30 = 0;
    let mut minlen: i64 = i64::MAX;
    let mut maxlen = 0;
    let mut len_vector: Vec<i64> = Vec::new();

    while !record.is_empty() {
        //let seq = record.seq();
        let len = record.seq().len() as i64; // here have to accomodate bigger numbers, as bases can get > 2^32
        
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
    let quart1 = quartiles(&mut len_vector, 1);
    let quart2 = quartiles(&mut len_vector, 2);
    let quart3 = quartiles(&mut len_vector, 3);
    let n50 = n50(&mut len_vector, 0.5); // use 0.1 for N90!!!
    let q20 = qual20 as f64 / bases as f64 * 100.0;
    let q30 = qual30 as f64 / bases as f64 * 100.0;

    println!("file\treads\tbases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent");
    println!("{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}", infile, reads, bases, minlen, maxlen, mean_len, quart1, quart2, quart3, n50, q20, q30);
}
