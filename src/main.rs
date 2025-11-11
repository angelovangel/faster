use bio::seq_analysis::gc::gc_content;
use modules::write_fastq;
use regex::{bytes::RegexSet, Regex};
use std::{fs, io::BufRead, io::BufReader, process, time::Duration};
use indicatif::{HumanCount, ProgressBar};
use kseq::parse_path;

//extern crate clap;
use clap::{App, Arg, ArgGroup};
// own functions
mod modules;

fn main() {
    let matches = App::new("faster")
                        .version("0.2.4")
                        .author("Angel Angelov <aangeloo@gmail.com>")
                        .about("fast statistics and more for 1 fastq file")
                        .arg(Arg::with_name("skip_header")
                            .short('s')
                            .long("skip_header")
                            .help("skip header in table output"))
                        .arg(Arg::with_name("table")
                            .short('t')
                            .long("table")
                            .help("Output a (tab separated) table with statistics"))
                        .arg(Arg::with_name("len")
                            .short('l')
                            .long("len")
                            .help("Output read lengths, one line per read"))
                        .arg(Arg::with_name("gc")
                            .short('g')
                            .long("gc")
                            .help("Output GC-content, one line per read"))
                        .arg(Arg::with_name("qscore")
                            .short('q')
                            .long("qscore")
                            .help("Output 'mean' read phred scores, one line per read. For this, the mean of the base probabilities for each read is calculated, and the result is converted back to a phred score"))
                        .arg(Arg::with_name("nx")
                            .long("nx")
                            .short('x')
                            .takes_value(true)
                            .help("Output NX value, provide the desired NX value as 0.5 for e.g. N50 [numeric]"))
                        .arg(Arg::with_name("qyield")
                            .long("qyield")
                            .short('y')
                            .takes_value(true)
                            .help("Percent bases with Q score of x or higher (use range 8..60)"))
                        .arg(Arg::with_name("sample")
                        //.short("s")
                            .long("sample")
                            .short('p')
                            .takes_value(true)
                            .help("Sub-sample sequences by proportion (0.0 to 1.0). Slow on large files!"))
                        .arg(Arg::with_name("filterl")
                            .short('f')
                            .long("filterl")
                            .takes_value(true)
                            .allow_hyphen_values(true) //important to parse negative integers
                            .help("Filter reads based on length - use positive integer to filter for reads LONGER than [integer] and negative integer to filter for reads that are SHORTER than [integer]"))
                        .arg(Arg::with_name("filterq")
                            .short('w')
                            .long("filterq")
                            .takes_value(true)
                            .allow_hyphen_values(true) //important to parse negative integers
                            .help("Filter reads based on 'mean' read quality - use positive integer to filter for reads with BETTER quality than [integer] and negative integer to filter for reads with WORSE wuality than [integer]. Use range 8..60 for qscores"))
                        .arg(Arg::with_name("trim_front")
                            .long("trim_front")
                            .short('a')
                            .takes_value(true)
                            .help("Trim all reads [integer] bases from the beginning"))
                        .arg(Arg::with_name("trim_tail")
                            .long("trim_tail")
                            .short('b')
                            .takes_value(true)
                            .help("Trim all reads [integer] bases from the end"))
                        .arg(Arg::with_name("regex_string")
                            .long("regex_string")
                            .short('r')
                            .takes_value(true)
                            .help("Output only reads whose id field matches a regex [string] pattern. See https://docs.rs/regex/1.4.2/regex/#functions"))
                        .arg(Arg::with_name("regex_file")
                            .long("regex_file")
                            .short('z')
                            .takes_value(true)
                            .help("Output only reads whose id field matches a regex [string] pattern. The regex patterns are read from a file, one line per pattern."))
                        .arg(Arg::with_name("INPUT")
                            .help("Paths to a fastq files, glob patterns and stdin can also be used")
                            .required(true)
                            .multiple_values(true)
                            .min_values(1)
                            .index(1))

                        // this group makes one and only one arg from the set required, avoid defining conflicts_with
                        .group(ArgGroup::with_name("group")
                        .required(true).args(&["table", "len", "gc", "qscore", "filterl", "filterq", "sample", "trim_front", "trim_tail", "regex_string", "regex_file", "nx", "qyield"]))
                        .get_matches();

    let infiles: Vec<&str> = matches.values_of("INPUT").unwrap().collect();
    
    // Check if the header for the table output needs to be printed once before the loop
    if matches.is_present("table") && !matches.is_present("skip_header") {
        println!("file\treads\tbases\tn_bases\tmin_len\tmax_len\tmean_len\tQ1\tQ2\tQ3\tN50\tQ20_percent\tQ30_percent");
    }


    for infile in infiles {
        let mut records = parse_path(&infile).unwrap();
    
        // case len
        if matches.is_present("len") {
            while let Some(record) = records.iter_record().unwrap() {
                println!("{}", record.len());
            }
            continue; // Go to the next file
            
        // case gc
        } else if matches.is_present("gc") {
            while let Some(record) = records.iter_record().unwrap() {
                let seq = record.seq().as_bytes();
                println!("{}", gc_content(seq));
            }
            continue; // Go to the next file

        // case qscore
        } else if matches.is_present("qscore") {
            while let Some(record) = records.iter_record().unwrap() {
                let mean_errorp = modules::qscore_probs(record.qual().as_bytes()) / record.seq().len() as f32;
                //println!("{:.4}", mean_errorp);
                println!("{:.4}", -10.0 * mean_errorp.log10());
            }
            continue; // Go to the next file

        // case filter length
        } else if matches.is_present("filterl") {
            let filterlen = matches
                .value_of("filterl")
                .unwrap()
                .trim()
                .parse::<i32>();
            // error on invalid input, rather than trying to guess
            match filterlen {
                Ok(x) => {
                    while let Some(record) = records.iter_record().unwrap() {
                        let seqlen = record.seq().len() as i32;
                        if x >= 0 {
                            if seqlen > x {
                                write_fastq(record);
                            }
                        } else if x < 0 {
                            if seqlen < x.abs() {
                                write_fastq(record);
                            }
                        }
                    }
                }
                Err(e) => eprintln!("Did you use an integer for filter? The error is: '{}'", e),
            }
            continue; // Go to the next file
        // case filter by qscore
        } else if matches.is_present("filterq"){
            let filterq = matches
                .value_of("filterq")
                .unwrap()
                .trim()
                .parse::<i32>()
                .expect("Failed to parse desired q value, please use integer between 8 and 60");
            match filterq {
                q if (-60..=60).contains(&q) => {
                    // do stuff
                    while let Some(record) = records.iter_record().unwrap() {
                        let mean_errorp = modules::qscore_probs(record.qual().as_bytes()) / record.seq().len() as f32;
                        let mean_qscore = -10.0 * mean_errorp.log10();
                        let q_f32 = q as f32;
                        if q >= 0 {
                            if mean_qscore > q_f32 {
                                write_fastq(record);
                            }
                        } else if q < 0 {
                            if mean_qscore < q_f32.abs() {
                                write_fastq(record);
                            }
                        }
                    }
                }
                _ => {
                    eprintln!("The q value should be between 10 and 60");
                    process::exit(1)
                }
            }
            continue; // Go to the next file
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

                        while let Some(record) = records.iter_record().unwrap() {
                            let len = record.seq().len() as i64;
                            lengths.push(len);                    
                        }
                        let nx = modules::get_nx(&mut lengths, 1.0 - nxvalue);
                        let nx_value = nxvalue * 100.0;
                        
                        println!("N{}\t{}", nx_value as i16, nx);
                    }
                    _ => {
                        eprintln!("The NX value should be between 0.1 and 1.0");
                        process::exit(0)
                    }
                }
            continue; // Go to the next file
        } else if matches.is_present("qyield") {
            let qvalue = matches
                .value_of("qyield")
                .unwrap()
                .trim()
                .parse::<u8>()
                .expect("Failed to parse desired q value");

                match qvalue {
                    x if (8..=60).contains(&x) => {
                        let mut bases: i64 = 0;
                        let mut qualx: i64 = 0;
                        // do work
                        while let Some(record) = records.iter_record().unwrap() {
                            let len = record.seq().len() as i64;
                            bases += len;
                            qualx += modules::get_qual_bases(record.qual().as_bytes(), 33 + qvalue); // 33 offset
                        }
                        let qx = qualx as f64 / bases as f64 * 100.0;
                        println!("Q{}\t{:.2}", qvalue, qx);
                    }
                    _ => {
                        eprintln!("The qyield value should be between 10 and 60");
                        process::exit(1)
                    }
                }
            continue; // Go to the next file
        
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
                        let nth = 1 as f32/fraction; // 1/fraction gives step_by
                        let mut recn = 0;
                        while let Some(record) = records.iter_record().unwrap() {
                            recn += 1;
                            if recn != nth as i32 {
                                continue;
                            }
                            recn = 0;
                            write_fastq(record);
                        }
                    }
                    _ => {
                        eprintln!("The subsample fraction should be between 0.0 and 1.0");
                        process::exit(0)
                    }
            }

            continue; // Go to the next file
            
        } else if matches.is_present("trim_front") {
            // parse trim value as usize
            let trimvalue = matches
                .value_of("trim_front")
                .unwrap()
                .trim()
                .parse::<usize>()
                .expect("failed to parse trim value!");

            while let Some(record) = records.iter_record().unwrap() {

                // the new sequence is trim..seq.len
                let newseq = &record.seq()[trimvalue..];
                let newqual = &record.qual()[trimvalue..];
                
                println!(
                    "{} {}\n{}\n{}\n{}", 
                    "@".to_string() + record.head(), record.des(), 
                    newseq, 
                    "+", 
                    newqual
                );
            }
            continue; // Go to the next file
        } else if matches.is_present("trim_tail") {
            let trimvalue = matches
                .value_of("trim_tail")
                .unwrap()
                .trim()
                .parse::<usize>()
                .expect("failed to parse trim value!");

            while let Some(record) = records.iter_record().unwrap() {
                
                // the new sequence is 0..seq.len - trim
                let trimright = record.len() - trimvalue;
                let newseq = &record.seq()[..trimright];
                let newqual = &record.qual()[..trimright];

                println!(
                    "{} {}\n{}\n{}\n{}", 
                    "@".to_string() + record.head(), record.des(), 
                    newseq, 
                    "+", 
                    newqual
                );
            }
            continue; // Go to the next file
        } else if matches.is_present("regex_string") {
            // parse string
            let string: &str = matches.value_of("regex_string").unwrap().trim();

            let re = Regex::new(string).expect("Failed to construct regex from string!");

            while let Some(record) = records.iter_record().unwrap() {
                let readid = record.head();
                if re.is_match(readid) {
                    write_fastq(record);  
                }
            }
            continue; // Go to the next file
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
            while let Some(record) = records.iter_record().unwrap() {
                let readid = record.head().as_bytes(); // as.bytes because RegexSet matches on bytes

                if re_set.is_match(readid) {
                    write_fastq(record);
                }
            }
            //println!("vector: {:?}", &revec);
            continue; // Go to the next file

        // case table (only runs if table is requested and none of the other single-task options were matched)
        } else if matches.is_present("table") {
            // normal case, output table
            let mut reads: i64 = 0;
            let mut bases: i64 = 0;
            let mut num_n: i32 = 0;
            let mut qual20: i64 = 0;
            let mut qual30: i64 = 0;
            let mut minlen: i64 = i64::MAX;
            let mut maxlen = 0;
            let mut len_vector: Vec<i64> = Vec::new();
            let pb = ProgressBar::new_spinner();
            pb.enable_steady_tick(Duration::from_millis(120));

            while let Some(record) = records.iter_record().unwrap() {
                //let seq = record.seq();
                let len = record.len() as i64; // here have to accomodate bigger numbers, as bases can get > 2^32

                reads += 1;
                bases += len;
                num_n += modules::get_n_bases(record.seq().as_bytes());
                qual20 += modules::get_qual_bases(record.qual().as_bytes(), 53); // 33 offset
                qual30 += modules::get_qual_bases(record.qual().as_bytes(), 63);
                minlen = len.min(minlen);
                maxlen = len.max(maxlen);
                len_vector.push(len);
                let message = format!("Processed reads: {}", HumanCount(reads as u64).to_string());
                pb.set_message(message);
            }

            let mean_len = modules::mean(&len_vector);
            let quart1 = modules::quartiles(&mut len_vector, 1);
            let quart2 = modules::quartiles(&mut len_vector, 2);
            let quart3 = modules::quartiles(&mut len_vector, 3);
            let n50 = modules::get_nx(&mut len_vector, 0.5); // use 0.1 for N90!!!
            let q20 = qual20 as f64 / bases as f64 * 100.0;
            let q30 = qual30 as f64 / bases as f64 * 100.0;
            pb.finish_and_clear();

            // The header is now printed once before the loop (see top of main)
            
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
            continue; // Go to the next file
        }
    }
}
// END