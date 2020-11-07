# fastq-stats

A (very) fast program for getting statistics about a fastq file, the way I need them, written in Rust

## Description


## Install

Compiled binaries are provided for x86_64 Linux, macOS and Windows - just download and run. If you need to run it on something else (your phone?!), you will have to compile it yourself (which is pretty easy though). Below is an example on how to setup a Rust toolchain and download and compile `fastq-stats`:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
git clone https://github.com/angelovangel/fastq-stats.git

cd fastq-stats
cargo build --release

# the binary is now under ./target/release/, run it like this:
./target/release/fastq-stats file.fastq.gz

```

## Usage

The program takes one fasq/fastq.gz file as an argument and outputs a tab-separated table with statistics to stdout.

```bash
fastq-stats file.fastq

# for many files, with parallel
parallel fastq-stats ::: /path/to/fastq/*.fastq.gz
```

## Performance

Not tested very thoroughly, but it should be 2-5x faster than `seqkit stats`, especially on gz files. This is not a fair test though, as seqkit offers many more functionalities. Below some benchmarks on a MacBook Pro with an 8-core 2.3 GHz Quad-Core Intel Core i5 and 8 GB RAM:

```bash
# for a Nanopore fastq file with 37k reads and 350M bases
fastq-stats
0.13s user 0.16s system 99% cpu 0.288 total

seqkit stats
0.30s user 0.20s system 99% cpu 0.504 total

# for an Illumina fastq.gz file 
fastq-stats
0.17s user 0.00s system 98% cpu 0.174 total

seqkit stats
1.05s user 0.04s system 101% cpu 1.075 total
```
