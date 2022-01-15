use super::{BwaReference, BwaSettings, PairedEndStats};
use std::ffi::{CStr, CString};
use std::ptr;
use std::sync::{Arc, Mutex};
use std::path::Path;
use rust_htslib::bam::record::Record;
use rust_htslib::bam::HeaderView;

/// A BWA aligner. Carries everything required to align
/// reads to a reference and generate BAM records.
pub struct BwaAligner {
    reference: BwaReference,
    header_view: Arc<Mutex<HeaderView>>,
    settings: BwaSettings,
    pe_stats: PairedEndStats,
}
// this is not automatically derived because of an interior
//   mutable pointer inside HeaderView. It _is_ mutated
//   by the Record::from_sam function, so guard it with a mutex
unsafe impl Sync for BwaAligner {}

impl BwaAligner {
    /// Load a BWA reference from the given path and use default BWA settings and paired-end structure.
    pub fn from_path<P: AsRef<Path>>(path: P) -> anyhow::Result<BwaAligner> {
        let bwa_ref = BwaReference::open(path)?;
        Ok(BwaAligner::new(
            bwa_ref,
            BwaSettings::new(),
            PairedEndStats::default(),
        ))
    }

    pub fn new(
        reference: BwaReference,
        settings: BwaSettings,
        pe_stats: PairedEndStats,
    ) -> BwaAligner {
        let header = reference.create_bam_header();
        let header_view = Arc::new(Mutex::new(HeaderView::from_header(&header)));

        BwaAligner {
            reference,
            header_view,
            settings,
            pe_stats,
        }
    }

    /// Align a read-pair to the reference.
    pub fn align_read_pair(
        &self,
        name: &[u8],
        r1: &[u8],
        q1: &[u8],
        r2: &[u8],
        q2: &[u8],
    ) -> (Vec<Record>, Vec<Record>) {
        let name = CString::new(name).unwrap();
        let raw_name = name.into_raw();

        // Prep input data -- need to make copy of reads since BWA will edit the strings in-place
        // FIXME - set an id -- used for a random hash
        let mut r1 = Vec::from(r1);
        let mut q1 = Vec::from(q1);
        let mut r2 = Vec::from(r2);
        let mut q2 = Vec::from(q2);

        let read1 = bwa_sys::bseq1_t {
            l_seq: r1.len() as i32,
            name: raw_name,
            seq: r1.as_mut_ptr() as *mut i8,
            qual: q1.as_mut_ptr() as *mut i8,
            comment: ptr::null_mut(),
            id: 0,
            sam: ptr::null_mut(),
        };

        let read2 = bwa_sys::bseq1_t {
            l_seq: r2.len() as i32,
            name: raw_name,
            seq: r2.as_mut_ptr() as *mut i8,
            qual: q2.as_mut_ptr() as *mut i8,
            comment: ptr::null_mut(),
            id: 0,
            sam: ptr::null_mut(),
        };

        let mut reads = [read1, read2];

        // Align the read pair. BWA will write the SAM data back to the bwa_sys::bseq1_t.sam field
        unsafe {
            let r = *(self.reference.bwt_data);
            let settings = self.settings.bwa_settings;
            bwa_sys::mem_process_seq_pe(
                &settings,
                r.bwt,
                r.bns,
                r.pac,
                reads.as_mut_ptr(),
                self.pe_stats.inner.as_ptr(),
            );
            let _ = CString::from_raw(raw_name);
        }

        // Parse the results from the SAM output & convert the htslib Records
        let sam1 = unsafe { CStr::from_ptr(reads[0].sam) };
        let sam2 = unsafe { CStr::from_ptr(reads[1].sam) };

        let recs1 = self.parse_sam_to_records(sam1.to_bytes());
        let recs2 = self.parse_sam_to_records(sam2.to_bytes());

        unsafe {
            libc::free(reads[0].sam as *mut libc::c_void);
            libc::free(reads[1].sam as *mut libc::c_void);
        }

        (recs1, recs2)
    }

    /// Align a single read to the reference.
    pub fn align_single_read(
        &self,
        name: &[u8],
        r1: &[u8],
        q1: &[u8]
    ) -> Vec<Record> {
        let name = CString::new(name).unwrap();
        let raw_name = name.into_raw();

        // Prep input data -- need to make copy of reads since BWA will edit the strings in-place
        // FIXME - set an id -- used for a random hash
        let mut r1 = Vec::from(r1);
        let mut q1 = Vec::from(q1);

        let read1 = bwa_sys::bseq1_t {
            l_seq: r1.len() as i32,
            name: raw_name,
            seq: r1.as_mut_ptr() as *mut i8,
            qual: q1.as_mut_ptr() as *mut i8,
            comment: ptr::null_mut(),
            id: 0,
            sam: ptr::null_mut(),
        };

        let mut reads = [read1];

        // Align the read pair. BWA will write the SAM data back to the bwa_sys::bseq1_t.sam field
        unsafe {
            let r = *(self.reference.bwt_data);
            let settings = self.settings.bwa_settings;
            bwa_sys::mem_process_seq_pe(
                &settings,
                r.bwt,
                r.bns,
                r.pac,
                reads.as_mut_ptr(),
                self.pe_stats.inner.as_ptr(),
            );
            let _ = CString::from_raw(raw_name);
        }

        // Parse the results from the SAM output & convert the htslib Records
        let sam1 = unsafe { CStr::from_ptr(reads[0].sam) };

        let recs1 = self.parse_sam_to_records(sam1.to_bytes());

        unsafe {
            libc::free(reads[0].sam as *mut libc::c_void);
        }

        recs1
    }

    fn parse_sam_to_records(&self, sam: &[u8]) -> Vec<Record> {
        let mut records = Vec::new();

        for slc in sam.split(|x| *x == b'\n') {
            if slc.len() > 0 {
                let record = {
                    let header_view = self.header_view.lock().unwrap();
                    Record::from_sam(&header_view, slc).unwrap()
                };
                records.push(record);
            }
        }

        records
    }
}
