use bio::io::{fastq, fastq::FastqRead};
use bio::seq_analysis::gc::gc_content;
use flate2::bufread;
//use rand::seq::IteratorRandom;
use regex::{bytes::RegexSet, Regex};
//use std::io::Read;
use std::{fs, io, io::BufRead, io::BufReader, process};

extern crate clap;
use clap::{App, Arg, ArgGroup};
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
    }
}

fn samplefq(path: &String, n: usize) {
    let records = fastq::Reader::new(get_fastq_reader(path))
        .records()
        .step_by(n);
    let mut writer = fastq::Writer::new(io::stdout());
    
    // sampled is a vector of result types
    //let sampled = records.choose_multiple(&mut rng, n);
    for record in records {
        // get the record
        let r = record.unwrap();
        
        writer
            .write_record(&r)
            .expect("Failed to write fastq record!");
    }
}

fn main() {
    let matches = App::new("faster")
                        .version("0.2.0")
                        .author("Angel Angelov <aangeloo@gmail.com>")
                        .about("fast statistics and more for 1 fastq file")

                        .arg(Arg::with_name("skip_header")
                            .short("s")
                            .long("skip_header")
                            .help("skip header in table output"))

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
                            .help("Output GC-content, one line per read"))

                        .arg(Arg::with_name("qscore")
                            .short("q")
                            .long("qscore")
                            .help("Output 'mean' read phred scores, one line per read. For this, the mean of the base probabilities for each read is calculated, and the result is converted back to a phred score"))

                        .arg(Arg::with_name("nx")
                            .long("nx")
                            .takes_value(true)
                            .help("Output NX value, provide the desired NX value as 0.5 for e.g. N50 [numeric]"))
                        .arg(Arg::with_name("qyield")
                            .long("qyield")
                            .takes_value(true)
                            .help("Percent bases with Q score of x or higher (use range 8..60)"))
                        .arg(Arg::with_name("sample")
                        //.short("s")
                            .long("sample")
                            .takes_value(true)
                            .help("Sub-sample sequences by proportion (0.0 to 1.0). Slow on large files!"))

                        .arg(Arg::with_name("filter")
                            //.short("f")
                            .long("filter")
                            .takes_value(true)
                            .allow_hyphen_values(true) //important to parse negative integers
                            .help("Filter reads based on length - use positive integer to filter for reads LONGER than [integer] and negative integer to filter for reads that are SHORTER than [integer]"))

                        .arg(Arg::with_name("trim_front")
                            .long("trim_front")
                            .takes_value(true)
                            .help("Trim all reads [integer] bases from the beginning"))

                        .arg(Arg::with_name("trim_tail")
                            .long("trim_tail")
                            .takes_value(true)
                            .help("Trim all reads [integer] bases from the end"))

                        .arg(Arg::with_name("regex_string")
                            .long("regex_string")
                            .takes_value(true)
                            .help("Output only reads whose description field matches a regex [string] pattern. See https://docs.rs/regex/1.4.2/regex/#functions"))
                        .arg(Arg::with_name("regex_file")
                            .long("regex_file")
                            .takes_value(true)
                            .help("Output only reads whose description field matches a regex [string] pattern. The regex patterns are read from a file, one line per pattern."))
                        .arg(Arg::with_name("INPUT")
                            .help("Path to a fastq file")
                            .required(true)
                            .index(1))

                        // this group makes one and only one arg from the set required, avoid defining conflicts_with
                        .group(ArgGroup::with_name("group")
                        .required(true).args(&["table", "len", "gc", "qscore", "filter", "sample", "trim_front", "trim_tail", "regex_string", "regex_file", "nx", "qyield"]))
                        .get_matches();
    //println!("Working on {}", matches.value_of("INPUT").unwrap());
    // read file
    // Here we can call .unwrap() because the argument is required.
    let infile = matches.value_of("INPUT").unwrap().to_string();

    let mut reader = fastq::Reader::new(get_fastq_reader(&infile));
    let mut record = fastq::Record::new();
    // let mut writer = fastq::Writer::new(io::stdout());
    // define writer in loop where necessary, otherwise undefined behaviour!

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

            //let qscore = modules::phred_gm( &record.qual() );
            //println!("{:.4}", qscore);
            
            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);

    // case filter
    } else if matches.is_present("filter") {
        let filterlen = matches.value_of("filter").unwrap().trim().parse::<i32>();
        //.unwrap_or(0);
        // error on invalid input, rather than trying to guess
        match filterlen {
            Ok(x) => {
                while !record.is_empty() {
                    let mut writer = fastq::Writer::new(io::stdout());
                    let seqlen = record.seq().len() as i32;
                    if x >= 0 {
                        if seqlen > x {
                            writer
                                .write_record(&record)
                                .expect("Failed to write fastq record!");
                        }
                    } else if x < 0 {
                        if seqlen < x.abs() {
                            writer
                                .write_record(&record)
                                .expect("Failed to write fastq record!")
                        }
                    }
                    reader
                        .read(&mut record)
                        .expect("Failed to parse fastq record!");
                }
            }
            Err(e) => eprintln!("Did you use an integer for filter? The error is: '{}'", e),
        }
        process::exit(0);

    // case nx    
    } else if matches.is_present("nx") {
        let nxvalue: f32 = matches
            .value_of("nx")
            .unwrap()
            .trim()
            .parse::<f32>()
            .expect("Failed to parse desired NX value");

            match nxvalue {
                x if (0.0..=1.0).contains(&x) => {
                    // do work
                    let mut lengths: Vec<i64> = Vec::new();

                    while !record.is_empty() {
                        let len = record.seq().len() as i64;
                        lengths.push(len);

                        reader
                            .read(&mut record)
                            .expect("Failed to parse fastq record!");
                    
                    }
                    let nx = modules::get_nx(&mut lengths, 1.0 - nxvalue);
                    let nx_value = nxvalue * 100.0;
                    
                    println!("N{}\t{}", nx_value as i16, nx);
                    process::exit(0)
                }
                _ => {
                    eprintln!("The NX value should be between 0.1 and 1.0");
                    process::exit(0)
                }
            }
    } else if matches.is_present("qyield") {
        let qvalue = matches
            .value_of("qyield")
            .unwrap()
            .trim()
            .parse::<u8>()
            .expect("Failed to parse desired QX value");

            match qvalue {
                x if (8..=60).contains(&x) => {
                    let mut bases: i64 = 0;
                    let mut qualx: i64 = 0;
                    // do work
                    while !record.is_empty() {
                        let len = record.seq().len() as i64;
                        bases += len;
                        qualx += modules::get_qual_bases(record.qual(), 33 + qvalue); // 33 offset

                        reader
                            .read(&mut record)
                            .expect("Failed to parse fastq record!")
                    }
                    let qx = qualx as f64 / bases as f64 * 100.0;
                    println!("Q{}\t{:.2}", qvalue, qx);
                    process::exit(0)
                }
                _ => {
                    eprintln!("The qyield value should be between 10 and 60");
                    process::exit(0)
                }
            }
    
    } else if matches.is_present("sample") {
        // parse fraction
        let fraction = matches
            .value_of("sample")
            .unwrap()
            .trim()
            .parse::<f32>()
            .expect("Failed to parse sample fraction value!");

            match fraction {
                // see <https://stackoverflow.com/a/58434531/8040734>
                x if (0.0..=1.0).contains(&x) => {

                    samplefq(&infile, (1 as f32/fraction) as usize); // 1/fraction gives step_by 
                    process::exit(0);
                }
                _ => {
                    eprintln!("The subsample fraction should be between 0.0 and 1.0");
                    process::exit(0)
                }
        }

        
    } else if matches.is_present("trim_front") {
        // parse trim value as usize
        let trimvalue = matches
            .value_of("trim_front")
            .unwrap()
            .trim()
            .parse::<usize>()
            .expect("failed to parse trim value!");

        while !record.is_empty() {
            // new writer?
            let mut writer = fastq::Writer::new(io::stdout());
            let check = record.check();
            if check.is_err() {
                panic!("I got a rubbish record!")
            }
            //let seqlen = record.seq().len();
            let id = record.id();
            let desc = record.desc();
            // the new sequence is trim..seq.len
            let newseq = &record.seq()[trimvalue..];
            let newqual = &record.qual()[trimvalue..];
            let newrec = fastq::Record::with_attrs(id, desc, newseq, newqual);

            writer
                .write_record(&newrec)
                .expect("Failed to write fastq record!");

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    } else if matches.is_present("trim_tail") {
        let trimvalue = matches
            .value_of("trim_tail")
            .unwrap()
            .trim()
            .parse::<usize>()
            .expect("failed to parse trim value!");

        while !record.is_empty() {
            let mut writer = fastq::Writer::new(io::stdout());
            let seqlen = record.seq().len();
            let id = record.id();
            let desc = record.desc();
            // the new sequence is 0..seq.len - trim
            let trimright = seqlen - trimvalue;
            let newseq = &record.seq()[..trimright];
            let newqual = &record.qual()[..trimright];
            let newrec = fastq::Record::with_attrs(id, desc, newseq, newqual);

            writer
                .write_record(&newrec)
                .expect("Failed to write fastq record!");

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    } else if matches.is_present("regex_string") {
        // parse string
        let string: &str = matches.value_of("regex_string").unwrap().trim();

        let re = Regex::new(string).expect("Failed to construct regex from string!");

        while !record.is_empty() {
            let mut writer = fastq::Writer::new(io::stdout());
            let desc = record.desc().unwrap();
            if re.is_match(desc) {
                writer
                    .write_record(&mut record)
                    .expect("Error writing fastq record!");
            }

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        process::exit(0);
    } else if matches.is_present("regex_file") {
        //parse file
        let refilepath = matches.value_of("regex_file").unwrap();
        let refile = fs::File::open(refilepath).expect("File not found!");
        let re_reader = BufReader::new(refile);

        // collect regex lines in a vec
        let mut revec = Vec::new();
        for line in re_reader.lines().map(|l| l.unwrap()) {
            //println!("line is: {}", line);
            revec.push(line);
        }

        let re_set = RegexSet::new(&revec).unwrap();
        // write record to stdout in case of match
        while !record.is_empty() {
            let mut writer = fastq::Writer::new(io::stdout());
            let desc = record.desc().unwrap().as_bytes(); // as.bytes because RegexSet matches on bytes

            if re_set.is_match(desc) {
                writer
                    .write_record(&mut record)
                    .expect("Error writing fastq record!");
            }

            reader
                .read(&mut record)
                .expect("Failed to parse fastq record!");
        }
        //println!("vector: {:?}", &revec);
        process::exit(0);
    }

    // normal case, output table
    let mut reads: i64 = 0;
    let mut bases: i64 = 0;
    let mut num_n: i32 = 0;
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
        num_n += modules::get_n_bases(record.seq());
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

    if !matches.is_present("skip_header") {
    println!("file\treads\tbases\tn_bases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent");
    }
    println!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{:.2}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}",
        infile,
        reads,
        bases,
        num_n,
        minlen,
        maxlen,
        mean_len,
        quart1,
        quart2,
        quart3,
        n50,
        q20,
        q30
    );
}
// END
