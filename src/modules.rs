// simple helper functions for calcuting mean, quartiles etc
use rayon::prelude::*;
use kseq::record::Fastx;
pub fn mean(numbers: &[i64]) -> f64 {
    numbers.par_iter().sum::<i64>() as f64 / numbers.len() as f64
}

pub fn quartiles(numbers: &mut [i64], q: i8) -> i64 {
    numbers.sort_unstable();

    match q {
        1 => {
            let index = numbers.len() / 4;
            return numbers[index];
        }
        2 => {
            let index = numbers.len() / 2;
            return numbers[index];
        }
        3 => {
            // avoid having to use f64
            let index1 = numbers.len() / 4;
            let index2 = numbers.len() / 2;
            return numbers[index1 + index2];
        }
        _ => 42, //:)
    }
    // first quartile
}

pub fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {
    numbers.sort_unstable();

    // half of the bases
    let halfsum = numbers.par_iter().sum::<i64>() as f32 * fraction; // f32 * f32

    // cumsum of the sorted vector
    let cumsum = numbers
        .iter()
        .scan(0, |sum, i| {
            *sum += i;
            Some(*sum)
        })
        .collect::<Vec<_>>();
    let n50_index = cumsum.par_iter().position_first(|&x| x > halfsum as i64).unwrap();

    numbers[n50_index]
}

// get number of bases with q >= value
pub fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
    let mut n = 0;
    for item in q {
        if *item >= qx {
            n += 1
        }
    }
    n
}

// get number of N bases
pub fn get_n_bases(seq: &[u8]) -> i32 {
    let mut n = 0;
    for s in seq {
        if s == &78u8 || s == &110u8 {
            n += 1;
        }
    }
    n
}

// to get mean of q scores from a record - first convert to prob, calc mean, then back to phred
// this fn reads phred and converts to probs and returns their sum
//
// see how seqkit is doing it
// https://github.com/shenwei356/bio/blob/1886d4a9eab7315f6f445595acbdc7bc3edf0e08/seq/seq.go#L727

pub fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum = 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred / 10.0);
        qprob_sum += prob
    }
    qprob_sum
}

pub fn write_fastq(rec: Fastx<'_>) {
    println!(
        "{} {}\n{}\n{}\n{}", 
        "@".to_string() + rec.head(), rec.des(), 
        rec.seq(), 
        "+", 
        rec.qual()
    );
}
// get geometric mean from phred scores
// pub fn phred_gm(q: &[u8]) -> f64 {
//     let mut phred_product = 0.;
    
//     for &item in q {                       
//     let phred = *&item as f64 - 33.0;

//             phred_product += phred.ln();

//     }
//     let result = phred_product / q.len() as f64;
//     result.exp()
// }

/// Map a DNA triplet to an index 0..63; returns None for ambiguous bases.
fn triplet_index(tri: &[u8]) -> Option<usize> {
    let encode = |b: u8| match b {
        b'A' | b'a' => Some(0usize),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    };
    let a = encode(tri[0])?;
    let b = encode(tri[1])?;
    let c = encode(tri[2])?;
    Some(a * 16 + b * 4 + c)
}

/// Calculate the DUST score for a sequence.
/// Uses windows of `window` bp (default 64). Returns the maximum window score
/// across the whole read. A score > 2.0 indicates low-complexity sequence.
pub fn dust_score(seq: &[u8], window: usize) -> f64 {
    if seq.len() < 3 {
        return 0.0;
    }
    let mut max_score: f64 = 0.0;
    for chunk in seq.chunks(window) {
        if chunk.len() < 3 {
            continue;
        }
        let mut counts = [0u32; 64]; // 4^3 possible DNA triplets
        let mut n_triplets: f64 = 0.0;
        for tri in chunk.windows(3) {
            if let Some(idx) = triplet_index(tri) {
                counts[idx] += 1;
                n_triplets += 1.0;
            }
        }
        if n_triplets < 2.0 {
            continue;
        }
        let score: f64 = counts
            .iter()
            .filter(|&&c| c > 0)
            .map(|&c| (c as f64 * (c as f64 - 1.0)) / 2.0)
            .sum::<f64>()
            / (n_triplets - 1.0);
        if score > max_score {
            max_score = score;
        }
    }
    max_score
}
