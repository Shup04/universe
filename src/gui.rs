use crate::simulation::{Universe, View};
use minifb::{Key, MouseMode, Window, WindowOptions};
use rand::Rng;

const SIM_WIDTH: usize = 400;
const SIM_HEIGHT: usize = 400;
const SCALE: usize = 1;
const WINDOW_WIDTH: usize = SIM_WIDTH * SCALE;
const WINDOW_HEIGHT: usize = SIM_HEIGHT * SCALE;

pub fn run() -> Result<(), minifb::Error> {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut stats_frame: u32 = 0;
    let mut last_title = String::from("Starting...");

    let mut window = Window::new(
        "Now Wave | LMB + | RMB - | Shift+L push right | Shift+R push left | H heat | T warp | P rho | C charge | X core | F momentum | Space reset",
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        WindowOptions::default(),
    )?;

    window.set_target_fps(60);

    let mut last_m = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if window.is_key_down(Key::Space) {
            universe.reset();
        }

        stats_frame = stats_frame.wrapping_add(1);
        let refresh_stats = stats_frame % 10 == 0;

        let is_l = window.get_mouse_down(minifb::MouseButton::Left);
        let is_r = window.get_mouse_down(minifb::MouseButton::Right);
        let shift = window.is_key_down(Key::LeftShift) || window.is_key_down(Key::RightShift);

        if !last_m {
            if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Discard) {
                let cx = mx as isize / SCALE as isize;
                let cy = my as isize / SCALE as isize;

                if shift {
                    if is_l {
                        universe.push_nearby(cx, cy, 0.8, 0.0);
                    }
                    if is_r {
                        universe.push_nearby(cx, cy, -0.8, 0.0);
                    }
                } else {
                    if is_l {
                        universe.spawn_vortex(cx, cy, 1.0);
                    }
                    if is_r {
                        universe.spawn_vortex(cx, cy, -1.0);
                    }
                }
            }
        }

        last_m = is_l || is_r;

        let steps = if window.is_key_down(Key::T) { 100 } else { 1 };

        for _ in 0..steps {
            universe.step();

            if window.is_key_down(Key::H) {
                let mut rng = rand::thread_rng();
                for i in 0..universe.u.len() {
                    if universe.rho[i] > 0.1 {
                        // Shake the vector field directly
                        universe.u[i] += rng.gen_range(-0.15..0.15);
                        universe.v[i] += rng.gen_range(-0.15..0.15);
                    }
                }
                universe.charge = universe.recompute_charge_from_uv(&universe.u, &universe.v);
            }
        }

        if refresh_stats && !window.is_key_down(Key::T) {
            let mut peaks = Vec::new();
            let mut search_rho = universe.rho.clone();

            for _ in 0..2 {
                let mut max_val = 0.0;
                let mut max_idx = 0;

                for (i, &val) in search_rho.iter().enumerate() {
                    if val > max_val {
                        max_val = val;
                        max_idx = i;
                    }
                }

                if max_val > 0.4 {
                    let px = (max_idx % SIM_WIDTH) as isize;
                    let py = (max_idx / SIM_WIDTH) as isize;

                    let mut sx = 0.0;
                    let mut sy = 0.0;
                    let mut total = 0.0;

                    for dy in -30..30 {
                        for dx in -30..30 {
                            let idx = universe.get_index(px + dx, py + dy);
                            let v = search_rho[idx];
                            if v > 0.1 {
                                sx += (px + dx) as f32 * v;
                                sy += (py + dy) as f32 * v;
                                total += v;
                                search_rho[idx] = 0.0;
                            }
                        }
                    }

                    if total > 0.0 {
                        peaks.push((sx / total, sy / total));
                    }
                }
            }

            let total_charge: f32 = universe.charge.iter().sum();
            let total_momentum: f32 = universe
                .px
                .iter()
                .zip(universe.py.iter())
                .map(|(px, py)| (px * px + py * py).sqrt())
                .sum();

            let mut max_core = 0.0_f32;
            for i in (0..universe.u.len()).step_by(2) {
                let x = (i % SIM_WIDTH) as isize;
                let y = (i / SIM_WIDTH) as isize;
                let core = universe.local_core_strength(x, y);
                if core > max_core {
                    max_core = core;
                }
            }

            let mut max_pair_overlap = 0.0_f32;
            for y in (0..SIM_HEIGHT).step_by(4) {
                for x in (0..SIM_WIDTH).step_by(4) {
                    let overlap = universe.local_pair_overlap(x as isize, y as isize);
                    if overlap > max_pair_overlap {
                        max_pair_overlap = overlap;
                    }
                }
            }

            let mut best_pos = (0.0_f32, 0isize, 0isize);
            let mut best_neg = (0.0_f32, 0isize, 0isize);

            for i in (0..universe.u.len()).step_by(2) {
                let x = (i % SIM_WIDTH) as isize;
                let y = (i / SIM_WIDTH) as isize;

                let core = universe.local_core_strength(x, y);
                let q = universe.charge[i];

                if q > 0.0 && core > best_pos.0 {
                    best_pos = (core, x, y);
                }
                if q < 0.0 && core > best_neg.0 {
                    best_neg = (core, x, y);
                }
            }

            let core_sep = if best_pos.0 > 0.0 && best_neg.0 > 0.0 {
                let dx = best_pos.1 as f32 - best_neg.1 as f32;
                let dy = best_pos.2 as f32 - best_neg.2 as f32;
                (dx * dx + dy * dy).sqrt()
            } else {
                0.0
            };

            last_title = if peaks.len() == 2 {
                let d =
                    ((peaks[0].0 - peaks[1].0).powi(2) + (peaks[0].1 - peaks[1].1).powi(2)).sqrt();
                format!(
            "Dist: {:.3} | CoreSep: {:.3} | +Core: {:.3} | -Core: {:.3} | PairOverlap: {:.3} | Charge: {:.3} | Mom: {:.3}",
            d,
            core_sep,
            best_pos.0,
            best_neg.0,
            max_pair_overlap,
            total_charge,
            total_momentum
        )
            } else {
                format!(
            "CoreSep: {:.3} | +Core: {:.3} | -Core: {:.3} | PairOverlap: {:.3} | Charge: {:.3} | Mom: {:.3}",
            core_sep,
            best_pos.0,
            best_neg.0,
            max_pair_overlap,
            total_charge,
            total_momentum
        )
            };
        }

        window.set_title(&last_title);

        let view = if window.is_key_down(Key::X) {
            View::Core
        } else if window.is_key_down(Key::F) {
            View::Momentum
        } else if window.is_key_down(Key::C) {
            View::Charge
        } else if window.is_key_down(Key::P) {
            View::Potential
        } else {
            View::Matter
        };
        let buffer = universe.frame(view);
        window.update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)?;
    }
    Ok(())
}
