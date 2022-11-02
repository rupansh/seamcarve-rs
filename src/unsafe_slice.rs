use std::cell::UnsafeCell;

pub struct UnsafeSlice<'a> {
    slice: &'a [UnsafeCell<f32>],
}
unsafe impl<'a> Send for UnsafeSlice<'a> {}
unsafe impl<'a> Sync for UnsafeSlice<'a> {}

impl<'a> UnsafeSlice<'a> {
    pub fn new(slice: &'a mut [f32]) -> Self {
        let ptr = slice as *mut [f32] as *const [UnsafeCell<f32>];
        Self {
            slice: unsafe { &*ptr },
        }
    }

    /// SAFETY: UB if read while another thread writes to the same index without
    /// synchronization
    pub unsafe fn read(&self, i: usize) -> f32 {
        let ptr = self.slice[i].get();
        *ptr
    }

    /// SAFETY: It is UB if two threads write to the same index without
    /// synchronization.
    pub unsafe fn write(&self, i: usize, value: f32) {
        let ptr = self.slice[i].get();
        *ptr = value;
    }
}
