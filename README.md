# fastq-stats

A (very) fast program for getting statistics about a fastq file, the way I need them, written in Rust

## Description

I wrote this program to get fast and accurate statistics about a fastq file, formatted as a tab-separated table. The motivation behind it:

- many of the tools out there are just wrong when it comes to calculating mean phred scores (yes, just taking the mean phred score is wrong, it should be the geometric mean)
- one simple executable doing one thing well, no dependencies
- the output can be easily tweaked as desired, just basic programming skills required
- written with speed in mind, and can be easily run in parallel

## Install

Compiled binaries are provided for x86_64 Linux, macOS and Windows - just download from the releases section and run. You will have to make the file executable (`chmod a+x fastq-stats...`) and for MacOS, allow running external apps in your security settings. If you need to run it on something else (your phone?!), you will have to compile it yourself (which is pretty easy though). Below is an example on how to setup a Rust toolchain and compile `fastq-stats`:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/angelovangel/fastq-stats.git

cd fastq-stats
cargo build --release

# the binary is now under ./target/release/, run it like this:
./target/release/fastq-stats /path/to/fastq/file.fastq.gz

```

## Usage and tweaking the output

The program takes one fasq/fastq.gz file as an argument and outputs a tab-separated table with statistics to stdout.

```bash
fastq-stats /path/to/fastq/file.fastq

# for many files, with parallel
parallel fastq-stats ::: /path/to/fastq/*.fastq.gz
```

The output is a tab-separated table with the following columns:   
`file   reads   bases   minlen  maxlen  av_len  median_len  N50 Q20_percent Q30_percent`

## Performance

Not tested very thoroughly, but it should be 2-3x faster than `seqkit stats -a`, especially for fastq.gz files. Below some benchmarks on a MacBook Pro with an 8-core 2.3 GHz Quad-Core Intel Core i5 and 8 GB RAM:

```bash
# for a Nanopore fastq file with 37k reads and 350M bases
fastq-stats
0.13s user 0.16s system 99% cpu 0.288 total

seqkit stats
0.30s user 0.20s system 99% cpu 0.504 total

wc -l
0.39s user 0.10s system 99% cpu 0.489 total

# for an Illumina fastq.gz file with 128k reads and 13M bases
fastq-stats
0.17s user 0.00s system 98% cpu 0.174 total

seqkit stats
1.05s user 0.04s system 101% cpu 1.075 total
```
