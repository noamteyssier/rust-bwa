
/// BWA settings object. Currently only default settings are enabled
pub struct BwaSettings {
    pub bwa_settings: bwa_sys::mem_opt_t,
}

impl BwaSettings {
    /// Create a `BwaSettings` object with default BWA parameters
    pub fn new() -> BwaSettings {
        let ptr = unsafe { bwa_sys::mem_opt_init() };
        let bwa_settings = unsafe { *ptr };
        unsafe { libc::free(ptr as *mut libc::c_void) };
        BwaSettings { bwa_settings }
    }

    /// Set alignment scores
    pub fn set_scores(
        mut self,
        matchp: i32,
        mismatch: i32,
        gap_open: i32,
        gap_extend: i32,
    ) -> BwaSettings {
        self.bwa_settings.a = matchp;
        self.bwa_settings.b = mismatch;
        self.bwa_settings.o_del = gap_open;
        self.bwa_settings.o_ins = gap_open;
        self.bwa_settings.e_del = gap_extend;
        self.bwa_settings.e_ins = gap_extend;

        unsafe {
            bwa_sys::bwa_fill_scmat(matchp, mismatch, self.bwa_settings.mat.as_mut_ptr());
        }
        self
    }

    /// Set clipping score penalties
    pub fn set_clip_scores(mut self, clip5: i32, clip3: i32) -> BwaSettings {
        self.bwa_settings.pen_clip5 = clip5;
        self.bwa_settings.pen_clip3 = clip3;
        self
    }

    /// Set unpaired read penalty
    pub fn set_unpaired(mut self, unpaired: i32) -> BwaSettings {
        self.bwa_settings.pen_unpaired = unpaired;
        self
    }

    /// Mark shorter splits as secondary
    pub fn set_no_multi(mut self) -> BwaSettings {
        self.bwa_settings.flag |= 0x10; // MEM_F_NO_MULTI
        self
    }
}
