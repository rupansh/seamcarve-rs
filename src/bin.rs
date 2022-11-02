use std::time::Instant;

use image::io::Reader;
use seamcarve::seam_carve;

fn main() {
    let im = Reader::open("./test.png").unwrap().decode().unwrap();
    let rgb8 = im.to_rgba8();
    let ins = Instant::now();
    let luma = seam_carve(&rgb8, 1000);
    println!("total elapsed {}", ins.elapsed().as_secs());
    luma.save("res.png").unwrap();
}
