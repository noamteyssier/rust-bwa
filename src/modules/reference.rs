use std::path::Path;
use std::ffi::{CStr, CString};
use rust_htslib::bam::header::{Header, HeaderRecord};

#[derive(Debug, thiserror::Error)]
#[error("{0}")]
pub struct ReferenceError(String);

/// A BWA reference object to perform alignments to.
/// Must be loaded from a BWA index created with `bwa index`
pub struct BwaReference {
    pub bwt_data: *const bwa_sys::bwaidx_t,
    contig_names: Vec<String>,
    contig_lengths: Vec<usize>,
}
unsafe impl Sync for BwaReference {}

impl BwaReference {
    /// Load a BWA reference from disk. Pass the fasta filename of the
    /// original reference as `path`
    pub fn open<P: AsRef<Path>>(path: P) -> Result<BwaReference, ReferenceError> {
        let idx_file = CString::new(path.as_ref().to_str().unwrap()).unwrap();
        let idx = unsafe { bwa_sys::bwa_idx_load(idx_file.as_ptr(), 0x7 as i32) }; // FIXME -- use BWA_IDX_ALL

        if idx.is_null() {
            return Err(ReferenceError(format!(
                "couldn't load reference: {:?}",
                path.as_ref()
            )));
        }

        let mut contig_names = Vec::new();
        let mut contig_lengths = Vec::new();
        let num_contigs = unsafe { (*(*idx).bns).n_seqs };

        for i in 0..num_contigs as isize {
            unsafe {
                let name = CStr::from_ptr((*(*(*idx).bns).anns.offset(i)).name);
                let sz = (*(*(*idx).bns).anns.offset(i)).len;

                let name_string = name.to_owned().into_string().unwrap();
                contig_names.push(name_string);
                contig_lengths.push(sz as usize)
            }
        }

        Ok(BwaReference {
            bwt_data: idx,
            contig_names,
            contig_lengths,
        })
    }

    pub fn create_bam_header(&self) -> Header {
        let mut header = Header::new();
        self.populate_bam_header(&mut header);
        header
    }

    pub fn populate_bam_header(&self, header: &mut Header) {
        for (ref contig_name, &len) in self.contig_names.iter().zip(self.contig_lengths.iter()) {
            add_ref_to_bam_header(header, &contig_name, len);
        }
    }
}

impl Drop for BwaReference {
    fn drop(&mut self) {
        unsafe {
            bwa_sys::bwa_idx_destroy(self.bwt_data as *mut bwa_sys::bwaidx_t);
        }
    }
}

fn add_ref_to_bam_header(header: &mut Header, seq_name: &str, seq_len: usize) {
    let mut header_rec = HeaderRecord::new(b"SQ");
    header_rec.push_tag(b"SN", &seq_name);
    header_rec.push_tag(b"LN", &seq_len);
    header.push_record(&header_rec);
}
