use packed_simd::u16x32;

pub const RGB_Y: u16x32 = u16x32::new(
    3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0, 3, 1, 4, 0,
);

pub const SOBEL_X: [f32; 9] = [1., 2., 1., 0., 0., 0., -1., -2., -1.];

pub const SOBEL_Y: [f32; 9] = [1., 0., -1., 2., 0., -2., 1., 0., -1.];

pub const BLOCK_WIDTH: usize = 80;
pub const BLOCK_HEIGHT: usize = BLOCK_WIDTH / 2;
