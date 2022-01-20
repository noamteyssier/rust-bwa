/* automatically generated by rust-bindgen 0.55.1 */
#![allow(non_snake_case, non_camel_case_types)]

pub const BWA_IDX_BWT: u32 = 1;
pub const BWA_IDX_BNS: u32 = 2;
pub const BWA_IDX_PAC: u32 = 4;
pub const BWA_IDX_ALL: u32 = 7;
pub type size_t = ::std::os::raw::c_ulong;
pub type bwtint_t = u64;
#[repr(C)]
#[derive(Copy, Clone)]
pub struct bwt_t {
    pub primary: bwtint_t,
    pub L2: [bwtint_t; 5usize],
    pub seq_len: bwtint_t,
    pub bwt_size: bwtint_t,
    pub bwt: *mut u32,
    pub cnt_table: [u32; 256usize],
    pub sa_intv: ::std::os::raw::c_int,
    pub n_sa: bwtint_t,
    pub sa: *mut bwtint_t,
}
#[test]
fn bindgen_test_layout_bwt_t() {
    assert_eq!(
        ::std::mem::size_of::<bwt_t>(),
        1120usize,
        concat!("Size of: ", stringify!(bwt_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bwt_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bwt_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).primary as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(primary)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).L2 as *const _ as usize },
        8usize,
        concat!("Offset of field: ", stringify!(bwt_t), "::", stringify!(L2))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).seq_len as *const _ as usize },
        48usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(seq_len)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).bwt_size as *const _ as usize },
        56usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(bwt_size)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).bwt as *const _ as usize },
        64usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(bwt)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).cnt_table as *const _ as usize },
        72usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(cnt_table)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).sa_intv as *const _ as usize },
        1096usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(sa_intv)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).n_sa as *const _ as usize },
        1104usize,
        concat!(
            "Offset of field: ",
            stringify!(bwt_t),
            "::",
            stringify!(n_sa)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwt_t>())).sa as *const _ as usize },
        1112usize,
        concat!("Offset of field: ", stringify!(bwt_t), "::", stringify!(sa))
    );
}
pub type __off_t = ::std::os::raw::c_long;
pub type __off64_t = ::std::os::raw::c_long;
pub type FILE = _IO_FILE;
pub type _IO_lock_t = ::std::os::raw::c_void;
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct _IO_marker {
    pub _next: *mut _IO_marker,
    pub _sbuf: *mut _IO_FILE,
    pub _pos: ::std::os::raw::c_int,
}
#[test]
fn bindgen_test_layout__IO_marker() {
    assert_eq!(
        ::std::mem::size_of::<_IO_marker>(),
        24usize,
        concat!("Size of: ", stringify!(_IO_marker))
    );
    assert_eq!(
        ::std::mem::align_of::<_IO_marker>(),
        8usize,
        concat!("Alignment of ", stringify!(_IO_marker))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_marker>()))._next as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_marker),
            "::",
            stringify!(_next)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_marker>()))._sbuf as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_marker),
            "::",
            stringify!(_sbuf)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_marker>()))._pos as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_marker),
            "::",
            stringify!(_pos)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct _IO_FILE {
    pub _flags: ::std::os::raw::c_int,
    pub _IO_read_ptr: *mut ::std::os::raw::c_char,
    pub _IO_read_end: *mut ::std::os::raw::c_char,
    pub _IO_read_base: *mut ::std::os::raw::c_char,
    pub _IO_write_base: *mut ::std::os::raw::c_char,
    pub _IO_write_ptr: *mut ::std::os::raw::c_char,
    pub _IO_write_end: *mut ::std::os::raw::c_char,
    pub _IO_buf_base: *mut ::std::os::raw::c_char,
    pub _IO_buf_end: *mut ::std::os::raw::c_char,
    pub _IO_save_base: *mut ::std::os::raw::c_char,
    pub _IO_backup_base: *mut ::std::os::raw::c_char,
    pub _IO_save_end: *mut ::std::os::raw::c_char,
    pub _markers: *mut _IO_marker,
    pub _chain: *mut _IO_FILE,
    pub _fileno: ::std::os::raw::c_int,
    pub _flags2: ::std::os::raw::c_int,
    pub _old_offset: __off_t,
    pub _cur_column: ::std::os::raw::c_ushort,
    pub _vtable_offset: ::std::os::raw::c_schar,
    pub _shortbuf: [::std::os::raw::c_char; 1usize],
    pub _lock: *mut _IO_lock_t,
    pub _offset: __off64_t,
    pub __pad1: *mut ::std::os::raw::c_void,
    pub __pad2: *mut ::std::os::raw::c_void,
    pub __pad3: *mut ::std::os::raw::c_void,
    pub __pad4: *mut ::std::os::raw::c_void,
    pub __pad5: size_t,
    pub _mode: ::std::os::raw::c_int,
    pub _unused2: [::std::os::raw::c_char; 20usize],
}
#[test]
fn bindgen_test_layout__IO_FILE() {
    assert_eq!(
        ::std::mem::size_of::<_IO_FILE>(),
        216usize,
        concat!("Size of: ", stringify!(_IO_FILE))
    );
    assert_eq!(
        ::std::mem::align_of::<_IO_FILE>(),
        8usize,
        concat!("Alignment of ", stringify!(_IO_FILE))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._flags as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_flags)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_read_ptr as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_read_ptr)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_read_end as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_read_end)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_read_base as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_read_base)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_write_base as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_write_base)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_write_ptr as *const _ as usize },
        40usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_write_ptr)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_write_end as *const _ as usize },
        48usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_write_end)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_buf_base as *const _ as usize },
        56usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_buf_base)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_buf_end as *const _ as usize },
        64usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_buf_end)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_save_base as *const _ as usize },
        72usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_save_base)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_backup_base as *const _ as usize },
        80usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_backup_base)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._IO_save_end as *const _ as usize },
        88usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_IO_save_end)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._markers as *const _ as usize },
        96usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_markers)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._chain as *const _ as usize },
        104usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_chain)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._fileno as *const _ as usize },
        112usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_fileno)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._flags2 as *const _ as usize },
        116usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_flags2)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._old_offset as *const _ as usize },
        120usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_old_offset)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._cur_column as *const _ as usize },
        128usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_cur_column)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._vtable_offset as *const _ as usize },
        130usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_vtable_offset)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._shortbuf as *const _ as usize },
        131usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_shortbuf)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._lock as *const _ as usize },
        136usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_lock)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._offset as *const _ as usize },
        144usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_offset)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>())).__pad1 as *const _ as usize },
        152usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(__pad1)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>())).__pad2 as *const _ as usize },
        160usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(__pad2)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>())).__pad3 as *const _ as usize },
        168usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(__pad3)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>())).__pad4 as *const _ as usize },
        176usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(__pad4)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>())).__pad5 as *const _ as usize },
        184usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(__pad5)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._mode as *const _ as usize },
        192usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_mode)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<_IO_FILE>()))._unused2 as *const _ as usize },
        196usize,
        concat!(
            "Offset of field: ",
            stringify!(_IO_FILE),
            "::",
            stringify!(_unused2)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bntann1_t {
    pub offset: i64,
    pub len: i32,
    pub n_ambs: i32,
    pub gi: u32,
    pub is_alt: i32,
    pub name: *mut ::std::os::raw::c_char,
    pub anno: *mut ::std::os::raw::c_char,
}
#[test]
fn bindgen_test_layout_bntann1_t() {
    assert_eq!(
        ::std::mem::size_of::<bntann1_t>(),
        40usize,
        concat!("Size of: ", stringify!(bntann1_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bntann1_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bntann1_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).offset as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(offset)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).len as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(len)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).n_ambs as *const _ as usize },
        12usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(n_ambs)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).gi as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(gi)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).is_alt as *const _ as usize },
        20usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(is_alt)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).name as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(name)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntann1_t>())).anno as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(bntann1_t),
            "::",
            stringify!(anno)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bntamb1_t {
    pub offset: i64,
    pub len: i32,
    pub amb: ::std::os::raw::c_char,
}
#[test]
fn bindgen_test_layout_bntamb1_t() {
    assert_eq!(
        ::std::mem::size_of::<bntamb1_t>(),
        16usize,
        concat!("Size of: ", stringify!(bntamb1_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bntamb1_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bntamb1_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntamb1_t>())).offset as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bntamb1_t),
            "::",
            stringify!(offset)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntamb1_t>())).len as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(bntamb1_t),
            "::",
            stringify!(len)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntamb1_t>())).amb as *const _ as usize },
        12usize,
        concat!(
            "Offset of field: ",
            stringify!(bntamb1_t),
            "::",
            stringify!(amb)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bntseq_t {
    pub l_pac: i64,
    pub n_seqs: i32,
    pub seed: u32,
    pub anns: *mut bntann1_t,
    pub n_holes: i32,
    pub ambs: *mut bntamb1_t,
    pub fp_pac: *mut FILE,
}
#[test]
fn bindgen_test_layout_bntseq_t() {
    assert_eq!(
        ::std::mem::size_of::<bntseq_t>(),
        48usize,
        concat!("Size of: ", stringify!(bntseq_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bntseq_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bntseq_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).l_pac as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(l_pac)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).n_seqs as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(n_seqs)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).seed as *const _ as usize },
        12usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(seed)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).anns as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(anns)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).n_holes as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(n_holes)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).ambs as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(ambs)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bntseq_t>())).fp_pac as *const _ as usize },
        40usize,
        concat!(
            "Offset of field: ",
            stringify!(bntseq_t),
            "::",
            stringify!(fp_pac)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bwaidx_t {
    pub bwt: *mut bwt_t,
    pub bns: *mut bntseq_t,
    pub pac: *mut u8,
    pub is_shm: ::std::os::raw::c_int,
    pub l_mem: i64,
    pub mem: *mut u8,
}
#[test]
fn bindgen_test_layout_bwaidx_t() {
    assert_eq!(
        ::std::mem::size_of::<bwaidx_t>(),
        48usize,
        concat!("Size of: ", stringify!(bwaidx_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bwaidx_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bwaidx_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).bwt as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(bwt)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).bns as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(bns)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).pac as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(pac)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).is_shm as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(is_shm)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).l_mem as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(l_mem)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bwaidx_t>())).mem as *const _ as usize },
        40usize,
        concat!(
            "Offset of field: ",
            stringify!(bwaidx_t),
            "::",
            stringify!(mem)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct bseq1_t {
    pub l_seq: ::std::os::raw::c_int,
    pub id: ::std::os::raw::c_int,
    pub name: *mut ::std::os::raw::c_char,
    pub comment: *mut ::std::os::raw::c_char,
    pub seq: *mut ::std::os::raw::c_char,
    pub qual: *mut ::std::os::raw::c_char,
    pub sam: *mut ::std::os::raw::c_char,
}
#[test]
fn bindgen_test_layout_bseq1_t() {
    assert_eq!(
        ::std::mem::size_of::<bseq1_t>(),
        48usize,
        concat!("Size of: ", stringify!(bseq1_t))
    );
    assert_eq!(
        ::std::mem::align_of::<bseq1_t>(),
        8usize,
        concat!("Alignment of ", stringify!(bseq1_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).l_seq as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(l_seq)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).id as *const _ as usize },
        4usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(id)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).name as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(name)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).comment as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(comment)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).seq as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(seq)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).qual as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(qual)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<bseq1_t>())).sam as *const _ as usize },
        40usize,
        concat!(
            "Offset of field: ",
            stringify!(bseq1_t),
            "::",
            stringify!(sam)
        )
    );
}
extern "C" {
    pub fn bwa_fill_scmat(a: ::std::os::raw::c_int, b: ::std::os::raw::c_int, mat: *mut i8);
}
extern "C" {
    pub fn bwa_idx_load(
        hint: *const ::std::os::raw::c_char,
        which: ::std::os::raw::c_int,
    ) -> *mut bwaidx_t;
}
extern "C" {
    pub fn bwa_idx_destroy(idx: *mut bwaidx_t);
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct mem_opt_t {
    pub a: ::std::os::raw::c_int,
    pub b: ::std::os::raw::c_int,
    pub o_del: ::std::os::raw::c_int,
    pub e_del: ::std::os::raw::c_int,
    pub o_ins: ::std::os::raw::c_int,
    pub e_ins: ::std::os::raw::c_int,
    pub pen_unpaired: ::std::os::raw::c_int,
    pub pen_clip5: ::std::os::raw::c_int,
    pub pen_clip3: ::std::os::raw::c_int,
    pub w: ::std::os::raw::c_int,
    pub zdrop: ::std::os::raw::c_int,
    pub max_mem_intv: u64,
    pub T: ::std::os::raw::c_int,
    pub flag: ::std::os::raw::c_int,
    pub min_seed_len: ::std::os::raw::c_int,
    pub min_chain_weight: ::std::os::raw::c_int,
    pub max_chain_extend: ::std::os::raw::c_int,
    pub split_factor: f32,
    pub split_width: ::std::os::raw::c_int,
    pub max_occ: ::std::os::raw::c_int,
    pub max_chain_gap: ::std::os::raw::c_int,
    pub n_threads: ::std::os::raw::c_int,
    pub chunk_size: ::std::os::raw::c_int,
    pub mask_level: f32,
    pub drop_ratio: f32,
    pub XA_drop_ratio: f32,
    pub mask_level_redun: f32,
    pub mapQ_coef_len: f32,
    pub mapQ_coef_fac: ::std::os::raw::c_int,
    pub max_ins: ::std::os::raw::c_int,
    pub max_matesw: ::std::os::raw::c_int,
    pub max_XA_hits: ::std::os::raw::c_int,
    pub max_XA_hits_alt: ::std::os::raw::c_int,
    pub mat: [i8; 25usize],
}
#[test]
fn bindgen_test_layout_mem_opt_t() {
    assert_eq!(
        ::std::mem::size_of::<mem_opt_t>(),
        168usize,
        concat!("Size of: ", stringify!(mem_opt_t))
    );
    assert_eq!(
        ::std::mem::align_of::<mem_opt_t>(),
        8usize,
        concat!("Alignment of ", stringify!(mem_opt_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).a as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(a)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).b as *const _ as usize },
        4usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(b)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).o_del as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(o_del)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).e_del as *const _ as usize },
        12usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(e_del)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).o_ins as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(o_ins)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).e_ins as *const _ as usize },
        20usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(e_ins)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).pen_unpaired as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(pen_unpaired)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).pen_clip5 as *const _ as usize },
        28usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(pen_clip5)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).pen_clip3 as *const _ as usize },
        32usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(pen_clip3)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).w as *const _ as usize },
        36usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(w)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).zdrop as *const _ as usize },
        40usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(zdrop)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_mem_intv as *const _ as usize },
        48usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_mem_intv)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).T as *const _ as usize },
        56usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(T)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).flag as *const _ as usize },
        60usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(flag)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).min_seed_len as *const _ as usize },
        64usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(min_seed_len)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).min_chain_weight as *const _ as usize },
        68usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(min_chain_weight)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_chain_extend as *const _ as usize },
        72usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_chain_extend)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).split_factor as *const _ as usize },
        76usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(split_factor)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).split_width as *const _ as usize },
        80usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(split_width)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_occ as *const _ as usize },
        84usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_occ)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_chain_gap as *const _ as usize },
        88usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_chain_gap)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).n_threads as *const _ as usize },
        92usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(n_threads)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).chunk_size as *const _ as usize },
        96usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(chunk_size)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).mask_level as *const _ as usize },
        100usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(mask_level)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).drop_ratio as *const _ as usize },
        104usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(drop_ratio)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).XA_drop_ratio as *const _ as usize },
        108usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(XA_drop_ratio)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).mask_level_redun as *const _ as usize },
        112usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(mask_level_redun)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).mapQ_coef_len as *const _ as usize },
        116usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(mapQ_coef_len)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).mapQ_coef_fac as *const _ as usize },
        120usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(mapQ_coef_fac)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_ins as *const _ as usize },
        124usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_ins)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_matesw as *const _ as usize },
        128usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_matesw)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_XA_hits as *const _ as usize },
        132usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_XA_hits)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).max_XA_hits_alt as *const _ as usize },
        136usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(max_XA_hits_alt)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_opt_t>())).mat as *const _ as usize },
        140usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_opt_t),
            "::",
            stringify!(mat)
        )
    );
}
#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct mem_pestat_t {
    pub low: ::std::os::raw::c_int,
    pub high: ::std::os::raw::c_int,
    pub failed: ::std::os::raw::c_int,
    pub avg: f64,
    pub std: f64,
}
#[test]
fn bindgen_test_layout_mem_pestat_t() {
    assert_eq!(
        ::std::mem::size_of::<mem_pestat_t>(),
        32usize,
        concat!("Size of: ", stringify!(mem_pestat_t))
    );
    assert_eq!(
        ::std::mem::align_of::<mem_pestat_t>(),
        8usize,
        concat!("Alignment of ", stringify!(mem_pestat_t))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_pestat_t>())).low as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_pestat_t),
            "::",
            stringify!(low)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_pestat_t>())).high as *const _ as usize },
        4usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_pestat_t),
            "::",
            stringify!(high)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_pestat_t>())).failed as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_pestat_t),
            "::",
            stringify!(failed)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_pestat_t>())).avg as *const _ as usize },
        16usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_pestat_t),
            "::",
            stringify!(avg)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<mem_pestat_t>())).std as *const _ as usize },
        24usize,
        concat!(
            "Offset of field: ",
            stringify!(mem_pestat_t),
            "::",
            stringify!(std)
        )
    );
}
extern "C" {
    pub fn mem_opt_init() -> *mut mem_opt_t;
}
extern "C" {
    pub fn mem_process_seq_pe(
        opt: *const mem_opt_t,
        bwt: *const bwt_t,
        bns: *const bntseq_t,
        pac: *const u8,
        seqs: *mut bseq1_t,
        pes: *const mem_pestat_t,
    );
}
extern "C" {
    pub fn bwa_idx_build(
        fa: *const ::std::os::raw::c_char,
        prefix: *const ::std::os::raw::c_char,
        algo_type: ::std::os::raw::c_int,
        block_size: ::std::os::raw::c_int,
    ) -> ::std::os::raw::c_int;
}
