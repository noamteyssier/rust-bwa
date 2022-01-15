// Copyright (c) 2017 10X Genomics, Inc. All rights reserved.

//! Rust-bwa provides a simple API wrapper around the BWA aligner.
//! Pass read-pair information in, and get Rust-htslib BAM records
//! back.
//!
//! ```
//! use bwa::BwaAligner;
//!
//! let bwa = BwaAligner::from_path(&"tests/test_ref.fa").unwrap();
//!
//! let r1 = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGG";
//! let q1 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//! let r2 = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAAT";
//! let q2 = b"2222222222222222222222222222222222222222222222222222222222222222222222222222222";
//!
//! let (r1_alns, _r2_alns) = bwa.align_read_pair(b"read_name", r1, q1, r2, q2);
//! println!("r1 mapping -- tid: {}, pos: {}", r1_alns[0].tid(), r1_alns[0].pos());
//! ```

#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

extern crate libc;
extern crate rust_htslib;
extern crate anyhow;
extern crate thiserror;

mod modules;
pub use modules::{
    BwaAligner,
    BwaSettings,
    BwaReference,
    PairedEndStats};


#[cfg(test)]
mod tests {
    use rust_htslib::bam::Record;
    use super::*;

    fn load_aligner() -> BwaAligner {
        let aln = BwaAligner::from_path("tests/test_ref.fa");
        aln.unwrap()
    }

    #[test]
    fn test_load_aligner() {
        let _ = load_aligner();
    }

    fn read_simple() -> [&'static [u8]; 5] {
        let name: &[u8] = b"@chr_727436_727956_3:0:0_1:0:0_0/1";
        let r1  : &[u8] = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGGCTGGCGCGGCTGATTAATGACATTCCTCTTCCCGGTACAACGGGCGTTGAGCGCGAACTTTTTCGCGCACT";
        let q1  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2  : &[u8] = b"TGCTGCGTAGCAGATCGACCCAGGCATTCCCTAGCGTGCTCATGCTCTGGCTGGTAAACGCACGGATGAGGGCAAAAATCACCGCAATCCCGCTGGCGGCAGAAAGAAAGTTTTGCACCGTTAAGCCCGCCATCTGGCTGAAATAGCTCA";
        let q2  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn read_single() -> [&'static [u8]; 3] {
        let name: &[u8] = b"@chr_727436_727956_3:0:0_1:0:0_0/1";
        let r1  : &[u8] = b"GATGGCTGCGCAAGGGTTCTTACTGATCGCCACGTTTTTACTGGTGTTAATGGTGCTGGCGCGTCCTTTAGGCAGCGGGCTGGCGCGGCTGATTAATGACATTCCTCTTCCCGGTACAACGGGCGTTGAGCGCGAACTTTTTCGCGCACT";
        let q1  : &[u8] = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1]
    }

    fn read_split() -> [&'static [u8]; 5] {
        let name = b"@chr_1561275_1561756_1:0:0_2:0:0_5c/1";
        let r1 = b"GCATCGATAAGCAGGTCAAATTCTCCCGTCATTATCACCTCTGCTACTTAAATTTCCCGCTTTATAAGCCGATTACGGCCTGGCATTACCCTATCCATAATTTAGGTGGGATGCCCGGTGCGTGGTTGGCAGATCCGCTGTTCTTTATTT";
        let q1 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        let r2 = b"TCATCGACCCAGGTATCATCGCGACGGGTACGATTACTGGCGAAGGTGAGAATGTTTAAAATCCAGCCGCCGAGTTTTTCAGCAATGGTCACCCATGACCAACCGGTGAACAACGTGAGGGCCGCTGCCCAAACGCATAGCAGCGCAATA";
        let q2 = b"222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222";

        [name, r1, q1, r2, q2]
    }

    fn align_read(r: [&[u8]; 5]) -> (Vec<Record>, Vec<Record>) {
        let bwa = load_aligner();
        bwa.align_read_pair(r[0], r[1], r[2], r[3], r[4])
    }

    fn align_single_read(r: [&[u8]; 3]) -> Vec<Record> {
        let bwa = load_aligner();
        bwa.align_single_read(r[0], r[1], r[2])
    }

    #[test]
    fn simple_align() {
        let (r1, r2) = align_read(read_simple());
        assert_eq!(r1[0].pos(), 727806);
        assert_eq!(r2[0].pos(), 727435);
    }

    #[test]
    fn split_align() {
        let (r1, r2) = align_read(read_split());
        assert_eq!(r1.len(), 2);
        assert_eq!(r1[0].pos(), 931375);
        assert_eq!(r1[1].pos(), 932605);
        assert_eq!(r2[0].pos(), 932937);
    }

    #[test]
    fn single_align() {
        let r1 = align_single_read(read_single());
        assert_eq!(r1.len(), 1);
        assert_eq!(r1[0].pos(), 727806);
    }

    #[test]
    fn header() {
        let reference = BwaReference::open("tests/test_ref.fa").unwrap();
        let hdr = b"@SQ\tSN:PhiX\tLN:5386\n@SQ\tSN:chr\tLN:4639675";
        assert_eq!(
            reference.create_bam_header().to_bytes().as_slice(),
            &hdr[..]
        );
    }
}
