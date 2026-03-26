use minifb::{Key, MouseMode, Window, WindowOptions};
use rand::Rng;
use std::f32::consts::PI;

const SIM_WIDTH: usize = 600;
const SIM_HEIGHT: usize = 600;
const SCALE: usize = 1;
const WINDOW_WIDTH: usize = SIM_WIDTH * SCALE;
const WINDOW_HEIGHT: usize = SIM_HEIGHT * SCALE;

pub struct Universe {
    width: usize,
    height: usize,
    theta: Vec<f32>,
    theta_prev: Vec<f32>,
    rho: Vec<f32>,
    charge: Vec<f32>,
}

fn wrap_angle(a: f32) -> f32 {
    let mut x = a;
    while x > PI {
        x -= 2.0 * PI;
    }
    while x < -PI {
        x += 2.0 * PI;
    }
    x
}

impl Universe {
    pub fn new(width: usize, height: usize) -> Self {
        let size = width * height;
        Self {
            width,
            height,
            theta: vec![0.0; size],
            theta_prev: vec![0.0; size],
            rho: vec![0.01; size],
            charge: vec![0.0; size],
        }
    }

    fn get_index(&self, x: isize, y: isize) -> usize {
        let wx = x.rem_euclid(self.width as isize) as usize;
        let wy = y.rem_euclid(self.height as isize) as usize;
        wy * self.width + wx
    }

    fn local_charge_from_theta(&self, theta_field: &[f32], x: isize, y: isize) -> f32 {
        let c = theta_field[self.get_index(x, y)];
        let r = theta_field[self.get_index(x + 1, y)];
        let ur = theta_field[self.get_index(x + 1, y - 1)];
        let u = theta_field[self.get_index(x, y - 1)];

        let d1 = wrap_angle(r - c);
        let d2 = wrap_angle(ur - r);
        let d3 = wrap_angle(u - ur);
        let d4 = wrap_angle(c - u);

        (d1 + d2 + d3 + d4) / (2.0 * PI)
    }

    fn recompute_charge(&self, theta_field: &[f32]) -> Vec<f32> {
        let mut next_charge = vec![0.0; self.width * self.height];
        for y in 0..self.height as isize {
            for x in 0..self.width as isize {
                let idx = self.get_index(x, y);
                next_charge[idx] = self.local_charge_from_theta(theta_field, x, y);
            }
        }
        next_charge
    }

    fn spawn_vortex(&mut self, cx: isize, cy: isize, sign: f32, vx: f32, vy: f32) {
        let radius: isize = 20;
        let sigma = 20.0_f32;

        for dy in -radius..=radius {
            for dx in -radius..=radius {
                let r2 = (dx * dx + dy * dy) as f32;
                if r2 <= (radius * radius) as f32 {
                    let idx = self.get_index(cx + dx, cy + dy);

                    let envelope = (-r2 / (2.0 * sigma * sigma)).exp();

                    // Smooth vortex angle
                    let angle = sign * (dy as f32).atan2(dx as f32);

                    // Write, do not add
                    self.theta[idx] = angle * envelope;
                    self.theta_prev[idx] = angle * envelope;

                    // Small velocity kick if wanted
                    self.theta_prev[idx] -= (vx * dx as f32 + vy * dy as f32) * 0.01 * envelope;

                    self.rho[idx] = self.rho[idx].max(0.9 * envelope);
                }
            }
        }
    }

    fn spawn_pair(&mut self, cx: isize, cy: isize, separation: isize) {
        self.spawn_vortex(cx - separation, cy,  1.0, 0.0, 0.0);
        self.spawn_vortex(cx + separation, cy, -1.0, 0.0, 0.0);
    }

    pub fn step(&mut self) {
        let mut next_t = vec![0.0; self.theta.len()];
        let mut next_r = vec![0.0; self.rho.len()];

        let base_speed = 0.4;

        for y in 0..self.height as isize {
            for x in 0..self.width as isize {
                let idx = self.get_index(x, y);

                let up = self.theta[self.get_index(x, y - 1)];
                let dn = self.theta[self.get_index(x, y + 1)];
                let lt = self.theta[self.get_index(x - 1, y)];
                let rt = self.theta[self.get_index(x + 1, y)];

                let ct = self.theta[idx];
                let pt = self.theta_prev[idx];
                let cr = self.rho[idx];

                // Time slows in dense matter, but never fully freezes
                let local_speed = base_speed / (1.0 + cr * 2.0);

                // Phase wave update
                let lap = (up - ct).sin() + (dn - ct).sin() + (lt - ct).sin() + (rt - ct).sin();
                let mut nt = (2.0 * ct) - pt + (local_speed * lap);

                // Mild universal friction
                nt = (nt - ct) * 0.99 + ct;

                // Local knot detector from phase variance
                let diff_up = (up - ct).sin().powi(2);
                let diff_dn = (dn - ct).sin().powi(2);
                let diff_lt = (lt - ct).sin().powi(2);
                let diff_rt = (rt - ct).sin().powi(2);
                let phase_variance = diff_up + diff_dn + diff_lt + diff_rt;

                // Signed local winding strength
                let q = self.local_charge_from_theta(&self.theta, x, y);
                let knot_strength = q.abs().min(1.0);

                // Blend variance and winding strength
                let target_r = (0.35 * phase_variance + 0.65 * knot_strength).clamp(0.0, 1.0);

                // Gentle density relaxation
                let nr = (cr * 0.96) + (target_r * 0.04);

                next_t[idx] = nt;
                next_r[idx] = nr.max(0.01);
            }
        }

        self.theta_prev = self.theta.clone();
        self.theta = next_t;
        self.rho = next_r;
        self.charge = self.recompute_charge(&self.theta);
    }
}

fn theta_rho_to_color(theta: f32, rho: f32, potential_view: bool) -> u32 {
    let brightness = if potential_view {
        (rho * 5.0).clamp(0.05, 1.0)
    } else if rho < 0.05 {
        0.0
    } else {
        rho
    };

    let red = ((theta.sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;
    let grn = (((theta + 2.0).sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;
    let blu = (((theta + 4.0).sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;

    (255 << 24) | (red << 16) | (grn << 8) | blu
}

fn charge_to_color(charge: f32, rho: f32) -> u32 {
    let q = charge.clamp(-1.0, 1.0);
    let strength = (q.abs() * rho.max(0.15)).clamp(0.0, 1.0);
    let c = (255.0 * strength) as u32;

    if q > 0.02 {
        // positive charge = red
        (255 << 24) | (c << 16)
    } else if q < -0.02 {
        // negative charge = blue
        (255 << 24) | c
    } else {
        // near zero = dim grey
        let g = (40.0 * rho.clamp(0.0, 1.0)) as u32;
        (255 << 24) | (g << 16) | (g << 8) | g
    }
}

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer = vec![0u32; WINDOW_WIDTH * WINDOW_HEIGHT];

    let mut window = Window::new(
        "Now Wave | LMB: spawn | RMB: push | H: heat | T: warp | P: rho | C: charge | Space: reset",
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        WindowOptions::default(),
    )
    .unwrap();

    window.set_target_fps(60);

    let mut last_m = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if window.is_key_down(Key::Space) {
            universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
        }

        let is_l = window.get_mouse_down(minifb::MouseButton::Left);
        let is_r = window.get_mouse_down(minifb::MouseButton::Right);

        if !last_m {
            if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Discard) {
                let cx = mx as isize / SCALE as isize;
                let cy = my as isize / SCALE as isize;

                if is_l {
                    // One single matter vortex
                    universe.spawn_vortex(cx, cy, 1.0, 0.0, 0.0);
                    universe.charge = universe.recompute_charge(&universe.theta);
                }

                if is_r {
                    // Intentional separated matter/antimatter pair
                    //universe.spawn_pair(cx, cy, 45);
                    universe.spawn_vortex(cx, cy, -1.0, 0.0, 0.0);
                    universe.charge = universe.recompute_charge(&universe.theta);
                }
            }
        }

        last_m = is_l || is_r;

        let steps = if window.is_key_down(Key::T) { 100 } else { 1 };

        for _ in 0..steps {
            universe.step();

            if window.is_key_down(Key::H) {
                let mut rng = rand::thread_rng();
                for i in 0..universe.theta.len() {
                    if universe.rho[i] > 0.1 {
                        universe.theta[i] += rng.gen_range(-0.1..0.1);
                    }
                }
                universe.charge = universe.recompute_charge(&universe.theta);
            }
        }

        // Simple centroid HUD for two strongest blobs
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

        if peaks.len() == 2 {
            let d = ((peaks[0].0 - peaks[1].0).powi(2) + (peaks[0].1 - peaks[1].1).powi(2)).sqrt();
            window.set_title(&format!(
                "Dist: {:.4} | Total charge: {:.4}",
                d, total_charge
            ));
        } else {
            window.set_title(&format!(
                "Total charge: {:.4}",
                total_charge
            ));
        }

        let show_potential = window.is_key_down(Key::P);
        let show_charge = window.is_key_down(Key::C);

        for i in 0..buffer.len() {
            buffer[i] = if show_charge {
                charge_to_color(universe.charge[i], universe.rho[i])
            } else {
                theta_rho_to_color(universe.theta[i], universe.rho[i], show_potential)
            };
        }

        window
            .update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)
            .unwrap();
    }
}
