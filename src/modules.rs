// simple helper functions for calcuting mean, quartiles etc


pub fn mean(numbers: &[i64]) -> f64 {
    numbers.iter().sum::<i64>() as f64 / numbers.len() as f64
}

pub fn quartiles(numbers: &mut [i64], q: i8) -> i64 {
        
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

pub fn get_nx(numbers: &mut [i64], fraction: f32) -> i64 {

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
pub fn get_qual_bases(q: &[u8], qx: u8) -> i64 {
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
pub fn qscore_probs(q: &[u8]) -> f32 {
    let mut qprob_sum= 0.0;
    for &item in q.iter() {
        let phred = *&item as f32 - 33.0;
        let prob = 10.0_f32.powf(-phred/10.0);
        qprob_sum += prob
    }
    qprob_sum
}