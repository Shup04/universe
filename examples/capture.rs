//! Deterministic capture using exactly the same diagnostic pixel functions as the GUI.
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use universe::{Universe, View};

fn main() -> std::io::Result<()> {
    let directory = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "output/capture".into());
    fs::create_dir_all(&directory)?;
    let mut sim = Universe::new(160, 120);
    sim.spawn_vortex(80, 60, 1.0);
    for _ in 0..40 {
        sim.step();
        assert!(sim.is_finite());
    }
    for (name, view) in [
        ("matter", View::Matter),
        ("potential", View::Potential),
        ("charge", View::Charge),
        ("core", View::Core),
        ("momentum", View::Momentum),
    ] {
        let mut file = BufWriter::new(File::create(format!("{directory}/{name}.ppm"))?);
        writeln!(file, "P6\n160 120\n255")?;
        for pixel in sim.frame(view) {
            file.write_all(&[(pixel >> 16) as u8, (pixel >> 8) as u8, pixel as u8])?;
        }
        file.flush()?;
    }
    println!("Saved five deterministic 160x120 diagnostic frames after 40 steps to {directory}");
    Ok(())
}
