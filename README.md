# seamcarve-rs

Parallel seam carving in rust  
WIP

## Usage

Depends on [image](https://crates.io/crates/image)  
Requires Rust Nightly

```rust
use image::io::Reader;
use seamcarve::seam_carve;

fn main() {
	let img = Reader::open("im.png").unwrap().decode().unwrap();
	let rgb = img.to_rgb8();
	let resized_luma = seam_carve(&rgb8, 1000);
	resized_luma.save("im_resized.png").unwrap();
}
```

## Performance Improvements

There are several ways to improve performance right now:

- Faster Edge Detection
- Cached Edge Detection (Recalculating energy for necessary pixels only)
- Improving the Seam container for better cache locality (we're using raw vec right now)

## TODO

- Horizontal Seams
- Return RGB image
