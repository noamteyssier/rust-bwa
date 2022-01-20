use std::ffi::CString;
use std::os::raw::c_int;
use anyhow;

pub fn build_reference(filename: &str, prefix: &str) -> anyhow::Result<()>{
    let idx_file = CString::new(filename)?;
    let idx_prefix = CString::new(prefix)?;
    unsafe { bwa_sys::bwa_idx_build(idx_file.as_ptr(),  idx_prefix.as_ptr(), 1, -1); }
    Ok(())
}


#[test]
fn test_build() {
    let filename = "tests/wat.fa";
    build_reference(filename, filename).expect("error building reference");
}
