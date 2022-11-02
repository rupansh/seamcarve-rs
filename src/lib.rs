#![feature(core_intrinsics)]
#![feature(test)]
extern crate test;

mod consts;
mod unsafe_slice;

use std::time::{Duration, Instant};

use image::{ImageBuffer, Luma, RgbaImage};
use packed_simd::{f32x16, u16x32, u16x4, u8x32, Cast, FromCast};

use consts::*;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};
use unsafe_slice::UnsafeSlice;

fn to_gray(image: &RgbaImage) -> Vec<f32> {
    let mut out = Vec::with_capacity(image.as_raw().capacity() / 4);
    unsafe {
        out.set_len(image.as_raw().len() / 4);
    }

    image
        .chunks_exact(32)
        .map(u8x32::from_slice_unaligned)
        .enumerate()
        .for_each(|(c, chunk)| {
            let gray: u16x32 = (u16x32::from_cast(chunk) * RGB_Y) >> 3;
            let grayt = unsafe { std::mem::transmute::<u16x32, [u16x4; 8]>(gray.cast()) };
            for i in 0..8 {
                out[c * 8 + i] = grayt[i].wrapping_sum() as f32;
            }
        });

    return out;
}

fn conv(im: &[f32], w: usize, h: usize, ker: [f32; 9]) -> Vec<f32> {
    // The kernel's input positions relative to the current pixel.
    let taps: &[(isize, isize)] = &[
        (-1, -1),
        (0, -1),
        (1, -1),
        (-1, 0),
        (0, 0),
        (1, 0),
        (-1, 1),
        (0, 1),
        (1, 1),
    ];

    let mut out = im.to_vec();

    let sum = match ker.iter().fold(0.0, |s, &item| s + item) {
        x if x == 0.0 => 1.0,
        sum => sum,
    };

    for y in 1..h - 1 {
        for x in 1..w - 1 {
            let mut t = 0.0;

            // TODO: There is no need to recalculate the kernel for each pixel.
            // Only a subtract and addition is needed for pixels after the first
            // in each row.
            for (&k, &(a, b)) in ker.iter().zip(taps.iter()) {
                let x0 = x as isize + a;
                let y0 = y as isize + b;

                let p = im[y0 as usize * w + x0 as usize];

                t += p * k;
            }

            let t = t / sum;
            out[y * w + x] = t;
        }
    }

    out
}

fn sobel_autovec(im: &[f32], w: usize, h: usize) -> Vec<f32> {
    let ex = conv(im, w, h, SOBEL_X);
    let ey = conv(im, w, h, SOBEL_Y);

    let mut out = im.to_vec();

    ex.chunks_exact(16)
        .zip(ey.chunks_exact(16))
        .map(|(cx, cy)| {
            (
                f32x16::from_slice_unaligned(cx),
                f32x16::from_slice_unaligned(cy),
            )
        })
        .zip(out.chunks_exact_mut(16))
        .for_each(|((cx, cy), out)| {
            let res = cx.abs() + cy.abs();
            res.write_to_slice_unaligned(out);
        });

    out
}

// Based on https://github.com/MPLLang/mpl/blob/0a12db7906867e38be90faacd64b600b4a9be078/examples/src/seam-carve/SCI.sml#L52
fn triangular_seam_search<'a>(energy: Vec<f32>, w: usize, h: usize, out_seam: &'a mut [f32]) {
    let seam_enrg = |out: &UnsafeSlice<'a>, i: usize, j: isize| {
        if j < 0 || j >= w as isize {
            f32::INFINITY
        } else {
            unsafe { out.read(i * w + j as usize) }
        }
    };

    let set_se = |out: &UnsafeSlice<'a>, e: &[f32], i: usize, j: isize| unsafe {
        out.write(
            i * w + j as usize,
            if i == 0 {
                0.0
            } else {
                e[i * w + j as usize]
                    + seam_enrg(out, i - 1, j)
                        .min(seam_enrg(out, i - 1, j - 1))
                        .min(seam_enrg(out, i - 1, j + 1))
            },
        )
    };

    let n_b = 1 + (w - 1) / BLOCK_WIDTH;

    let triangle =
        |out: &UnsafeSlice<'a>, e: &[f32], i: usize, jm: usize, jm_f: fn(usize) -> isize| {
            for k in 0..(h - i).min(BLOCK_HEIGHT) {
                let lo = (jm as isize - jm_f(k)).max(0);
                let hi = (jm as isize + jm_f(k)).min(w as isize);
                for j in lo..hi {
                    set_se(out, e, i + k, j);
                }
            }
        };

    let utriang = |out, e, i, jm| triangle(out, e, i, jm, |k| BLOCK_HEIGHT as isize - k as isize);
    let ltriang = |out, e, i, jm| triangle(out, e, i, jm, |k| k as isize + 1);

    // SAFETY: Disjoint index access
    // see https://shwestrick.github.io/2020/07/29/seam-carve.html for an explanation
    let set_strip = |out, i| {
        (0..n_b)
            .into_par_iter()
            .for_each(|b| utriang(out, &energy, i, b * BLOCK_WIDTH + BLOCK_HEIGHT));
        (0..n_b)
            .into_par_iter()
            .for_each(|b| ltriang(out, &energy, i + 1, b * BLOCK_WIDTH));
    };

    let mut i = 0;
    let uout = UnsafeSlice::new(out_seam);
    while i < h {
        set_strip(&uout, i);
        i += BLOCK_HEIGHT + 1;
    }
}

fn remove_min_seam(img: &[f32], seams: &[f32], w: usize, h: usize) -> Vec<f32> {
    let mut nimg = vec![0.; img.len() - h];
    let idx_min_dob = |(j1, m1): (isize, f32), (j2, m2)| if m1 > m2 { (j2, m2) } else { (j1, m1) };
    let idx_min_trip = |a, b, c| idx_min_dob(a, idx_min_dob(b, c));
    let m_e = |i, j| {
        if j < 0 || j >= w as isize {
            f32::INFINITY
        } else {
            seams[i * w + j as usize]
        }
    };

    let (j_m, _) = (0..w as isize)
        .into_par_iter()
        .map(|j| (j, m_e(h - 1, j)))
        .reduce(|| (-1 as isize, f32::INFINITY), |r, j| idx_min_dob(r, j));

    let mut j = j_m;
    for i in (1..h).rev() {
        let id = i * w;
        let nid = i * (w - 1);
        let fid = id + j as usize;
        let fnid = nid + j as usize;
        (&mut nimg[nid..fnid]).copy_from_slice(&img[id..fid]);
        (&mut nimg[fnid..nid + w - 1]).copy_from_slice(&img[fid + 1..id + w]);

        (j, _) = idx_min_trip(
            (j, m_e(i - 1, j)),
            (j - 1, m_e(i - 1, j - 1)),
            (j + 1, m_e(i - 1, j + 1)),
        );
    }
    (&mut nimg[0..j as usize]).copy_from_slice(&img[0..j as usize]);
    (&mut nimg[j as usize..w - 1]).copy_from_slice(&img[j as usize + 1..w]);

    return nimg;
}

pub fn seam_carve(image: &RgbaImage, seams_to_remove: usize) -> ImageBuffer<Luma<u8>, Vec<u8>> {
    let h = image.height() as usize;
    let mut w = image.width() as usize;
    let mut luma = to_gray(&image);
    let mut elapsed = Duration::ZERO;
    // TODO: switch structure for better cache locality
    let mut seam_energies = vec![0.; w * h];
    for _ in 0..seams_to_remove {
        let ins = Instant::now();
        let sobel = sobel_autovec(&luma, w, h);
        elapsed += ins.elapsed();

        triangular_seam_search(sobel, w, h, &mut seam_energies);
        luma = remove_min_seam(&luma, &seam_energies, w, h);

        w -= 1;
    }

    println!("elapsed {}", elapsed.as_secs());

    let luma = ImageBuffer::from_vec(
        w as u32,
        h as u32,
        luma.into_iter().map(|v| v as u8).collect(),
    )
    .unwrap();

    return luma;
}
