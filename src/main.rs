use minifb::{Key, MouseMode, Window, WindowOptions};
use rand::Rng;
use std::f32::consts::PI;

const SIM_WIDTH: usize = 400;
const SIM_HEIGHT: usize = 400;
const SCALE: usize = 1;
const WINDOW_WIDTH: usize = SIM_WIDTH * SCALE;
const WINDOW_HEIGHT: usize = SIM_HEIGHT * SCALE;

pub struct Universe {
    width: usize, height: usize,
    rho: Vec<f32>,
    charge: Vec<f32>,
    px: Vec<f32>, py: Vec<f32>,
    u: Vec<f32>, u_prev: Vec<f32>,
    v: Vec<f32>, v_prev: Vec<f32>,
    w: Vec<f32>, w_prev: Vec<f32>,
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
            rho: vec![0.01; size],
            charge: vec![0.0; size],
            px: vec![0.0; size],
            py: vec![0.0; size],
            u: vec![0.0; size], u_prev: vec![0.0; size],
            v: vec![0.0; size], v_prev: vec![0.0; size],
            w: vec![-1.0; size], w_prev: vec![-1.0; size],
        }
    }

    fn get_index(&self, x: isize, y: isize) -> usize {
        let wx = x.rem_euclid(self.width as isize) as usize;
        let wy = y.rem_euclid(self.height as isize) as usize;
        wy * self.width + wx
    }

fn recompute_charge_from_uv(&self, u_field: &[f32], v_field: &[f32]) -> Vec<f32> {
        let mut next_charge = vec![0.0; self.width * self.height];

        for y in 0..self.height as isize {
            for x in 0..self.width as isize {
                let idx = self.get_index(x, y);

                // --- THE SQUELCH FILTER (Killing the Quantum Foam) ---
                // If U and V are microscopic, we are in the dead vacuum. 
                // Skip the atan2 calculation so we don't amplify the static.
                if u_field[idx].abs() < 0.01 && v_field[idx].abs() < 0.01 {
                    next_charge[idx] = 0.0;
                    continue;
                }

                let c_theta  = v_field[self.get_index(x, y)].atan2(u_field[self.get_index(x, y)]);
                let r_theta  = v_field[self.get_index(x + 1, y)].atan2(u_field[self.get_index(x + 1, y)]);
                let ur_theta = v_field[self.get_index(x + 1, y - 1)].atan2(u_field[self.get_index(x + 1, y - 1)]);
                let u_theta  = v_field[self.get_index(x, y - 1)].atan2(u_field[self.get_index(x, y - 1)]);

                let d1 = wrap_angle(r_theta - c_theta);
                let d2 = wrap_angle(ur_theta - r_theta);
                let d3 = wrap_angle(u_theta - ur_theta);
                let d4 = wrap_angle(c_theta - u_theta);

                next_charge[idx] = (d1 + d2 + d3 + d4) / (2.0 * std::f32::consts::PI);
            }
        }

        next_charge
    }

    fn amp_at(&self, idx: usize) -> f32 {
        // Now accurately measures the full 3D sphere
        (self.u[idx].powi(2) + self.v[idx].powi(2) + self.w[idx].powi(2)).sqrt()
    }

fn local_phase_variance(&self, x: isize, y: isize) -> f32 {
        let idx = self.get_index(x, y);

        let u_ct = self.u[idx]; let v_ct = self.v[idx]; let w_ct = self.w[idx];

        let u_up = self.u[self.get_index(x, y - 1)]; let u_dn = self.u[self.get_index(x, y + 1)];
        let u_lt = self.u[self.get_index(x - 1, y)]; let u_rt = self.u[self.get_index(x + 1, y)];

        let v_up = self.v[self.get_index(x, y - 1)]; let v_dn = self.v[self.get_index(x, y + 1)];
        let v_lt = self.v[self.get_index(x - 1, y)]; let v_rt = self.v[self.get_index(x + 1, y)];

        // Gather W for telemetry
        let w_up = self.w[self.get_index(x, y - 1)]; let w_dn = self.w[self.get_index(x, y + 1)];
        let w_lt = self.w[self.get_index(x - 1, y)]; let w_rt = self.w[self.get_index(x + 1, y)];

        let diff_up = (u_up - u_ct).powi(2) + (v_up - v_ct).powi(2) + (w_up - w_ct).powi(2);
        let diff_dn = (u_dn - u_ct).powi(2) + (v_dn - v_ct).powi(2) + (w_dn - w_ct).powi(2);
        let diff_lt = (u_lt - u_ct).powi(2) + (v_lt - v_ct).powi(2) + (w_lt - w_ct).powi(2);
        let diff_rt = (u_rt - u_ct).powi(2) + (v_rt - v_ct).powi(2) + (w_rt - w_ct).powi(2);

        diff_up + diff_dn + diff_lt + diff_rt
    }

fn local_knot_strength(&self, x: isize, y: isize) -> f32 {
    let u_up = self.u[self.get_index(x, y - 1)];
    let u_dn = self.u[self.get_index(x, y + 1)];
    let u_lt = self.u[self.get_index(x - 1, y)];
    let u_rt = self.u[self.get_index(x + 1, y)];

    let v_up = self.v[self.get_index(x, y - 1)];
    let v_dn = self.v[self.get_index(x, y + 1)];
    let v_lt = self.v[self.get_index(x - 1, y)];
    let v_rt = self.v[self.get_index(x + 1, y)];

    let du_dx = (u_rt - u_lt) * 0.5;
    let dv_dy = (v_dn - v_up) * 0.5;
    let du_dy = (u_dn - u_up) * 0.5;
    let dv_dx = (v_rt - v_lt) * 0.5;

    let jacobian = du_dx * dv_dy - du_dy * dv_dx;
    (jacobian * 5.0).abs().min(1.0)
}

fn local_core_strength(&self, x: isize, y: isize) -> f32 {
    let idx = self.get_index(x, y);
    let amp = self.amp_at(idx);
    let knot = self.local_knot_strength(x, y);
    let charge = self.charge[idx].abs().min(1.0);

    // core = strong topology + strong charge + low amplitude center
    (0.45 * knot + 0.45 * charge + 0.10 * (1.0 - amp)).clamp(0.0, 1.0)
}

fn signed_core_strength(&self, x: isize, y: isize) -> f32 {
    self.local_core_strength(x, y) * self.charge[self.get_index(x, y)]
}

fn local_pos_core(&self, x: isize, y: isize) -> f32 {
    let idx = self.get_index(x, y);
    let q = self.charge[idx].max(0.0);
    self.local_core_strength(x, y) * q.min(1.0)
}

fn local_neg_core(&self, x: isize, y: isize) -> f32 {
    let idx = self.get_index(x, y);
    let q = (-self.charge[idx]).max(0.0);
    self.local_core_strength(x, y) * q.min(1.0)
}

fn local_pair_overlap(&self, x: isize, y: isize) -> f32 {
    let mut pos_sum = 0.0;
    let mut neg_sum = 0.0;

    for dy in -5..=5 {
        for dx in -5..=5 {
            let core = self.local_core_strength(x + dx, y + dy);
            let q = self.charge[self.get_index(x + dx, y + dy)];

            if q > 0.0 {
                pos_sum += core * q.min(1.0);
            } else if q < 0.0 {
                neg_sum += core * (-q).min(1.0);
            }
        }
    }

    (0.01 * pos_sum * neg_sum).min(1.0)
}

fn local_overlap_pressure(&self, x: isize, y: isize) -> f32 {
    let idx = self.get_index(x, y);
    let amp = self.amp_at(idx);

    let phase = (self.local_phase_variance(x, y) * 0.25).min(1.0);

    // rises when the oscillator collapses or gets highly compressed
    (0.7 * (1.0 - amp).max(0.0) + 0.3 * phase).clamp(0.0, 1.0)
}

fn spawn_vortex(&mut self, cx: isize, cy: isize, sign: f32) {
    let w = self.width as isize;
    let h = self.height as isize;
    
    let lambda = 12.0_f32;
    let lambda_sq = lambda * lambda;

    // --- THE QUARANTINE ZONE ---
    let max_r = 160.0_f32; 
    let fade_start = 80.0_f32;

    for y in 0..h {
        for x in 0..w {
            let idx = self.get_index(x, y);

            let mut dx = x - cx;
            let mut dy = y - cy;
            if dx > w/2 { dx -= w; } else if dx < -w/2 { dx += w; }
            if dy > h/2 { dy -= h; } else if dy < -h/2 { dy += h; }
            
            let fx = dx as f32;
            let fy = dy as f32;
            let r_sq = fx * fx + fy * fy;
            let r = r_sq.sqrt();

            let theta = sign * fy.atan2(fx);

            // 1. Calculate the exact Belavin-Polyakov Instanton
            let denom = lambda_sq + r_sq;
            let knot_w = (lambda_sq - r_sq) / denom;
            let uv_mag = (2.0 * lambda * r) / denom;
            let knot_u = uv_mag * theta.cos();
            let knot_v = uv_mag * theta.sin();

            // 2. The Isolation Envelope (Flattening the Crease)
            let mut envelope = 1.0_f32;
            if r > max_r {
                envelope = 0.0;
            } else if r > fade_start {
                let t = (r - fade_start) / (max_r - fade_start);
                // The Smoothstep Polynomial
                envelope = 1.0 - (t * t * (3.0 - 2.0 * t)); 
            }

            // 3. Blend smoothly with the dead vacuum (0, 0, -1)
            let blended_u = knot_u * envelope;
            let blended_v = knot_v * envelope;
            let blended_w = knot_w * envelope + (-1.0) * (1.0 - envelope);

            // 4. Normalize! (Forces amplitude back to exactly 1.0 to silence the Mexican Hat)
            let mag = (blended_u.powi(2) + blended_v.powi(2) + blended_w.powi(2)).sqrt();
            let final_u = blended_u / mag;
            let final_v = blended_v / mag;
            let final_w = blended_w / mag;

            // 5. Write the safe geometry to the grid
            self.u[idx] = final_u;
            self.v[idx] = final_v;
            self.w[idx] = final_w;
            self.u_prev[idx] = final_u;
            self.v_prev[idx] = final_v;
            self.w_prev[idx] = final_w;
        }
    }
}

    fn push_nearby(&mut self, cx: isize, cy: isize, push_x: f32, push_y: f32) {
        let radius: isize = 28;

        for dy in -radius..=radius {
            for dx in -radius..=radius {
                let r2 = (dx * dx + dy * dy) as f32;
                if r2 <= (radius * radius) as f32 {
                    let idx = self.get_index(cx + dx, cy + dy);

                    if self.rho[idx] > 0.15 {
                        let envelope = (-r2 / 180.0).exp();
                        self.px[idx] += push_x * envelope;
                        self.py[idx] += push_y * envelope;
                    }
                }
            }
        }
    }

pub fn step(&mut self) {
    let mut next_u = vec![0.0; self.u.len()];
    let mut next_u_prev = vec![0.0; self.u.len()];
    let mut next_v = vec![0.0; self.v.len()];
    let mut next_v_prev = vec![0.0; self.v.len()];
    let mut next_w = vec![0.0; self.w.len()];
    let mut next_w_prev = vec![0.0; self.w.len()];
    let mut next_r = vec![0.0; self.rho.len()];
    let mut next_px = vec![0.0; self.px.len()];
    let mut next_py = vec![0.0; self.py.len()];

    let base_speed = 0.15;

    for y in 0..self.height as isize {
        for x in 0..self.width as isize {
            let idx = self.get_index(x, y);

            // 1. Current and Previous State
            let u_ct = self.u[idx]; let u_pt = self.u_prev[idx];
            let v_ct = self.v[idx]; let v_pt = self.v_prev[idx];
            let w_ct = self.w[idx]; let w_pt = self.w_prev[idx];

            // 2. High-Precision Gradients
            let u_rt = self.u[self.get_index(x+1,y)]; let u_lt = self.u[self.get_index(x-1,y)];
            let u_dn = self.u[self.get_index(x,y+1)]; let u_up = self.u[self.get_index(x,y-1)];
            let u_ur = self.u[self.get_index(x+1,y-1)]; let u_ul = self.u[self.get_index(x-1,y-1)];
            let u_dr = self.u[self.get_index(x+1,y+1)]; let u_dl = self.u[self.get_index(x-1,y+1)];
            
            let v_rt = self.v[self.get_index(x+1,y)]; let v_lt = self.v[self.get_index(x-1,y)];
            let v_dn = self.v[self.get_index(x,y+1)]; let v_up = self.v[self.get_index(x,y-1)];
            let v_ur = self.v[self.get_index(x+1,y-1)]; let v_ul = self.v[self.get_index(x-1,y-1)];
            let v_dr = self.v[self.get_index(x+1,y+1)]; let v_dl = self.v[self.get_index(x-1,y+1)];

            let w_rt = self.w[self.get_index(x+1,y)]; let w_lt = self.w[self.get_index(x-1,y)];
            let w_dn = self.w[self.get_index(x,y+1)]; let w_up = self.w[self.get_index(x,y-1)];
            let w_ur = self.w[self.get_index(x+1,y-1)]; let w_ul = self.w[self.get_index(x-1,y-1)];
            let w_dr = self.w[self.get_index(x+1,y+1)]; let w_dl = self.w[self.get_index(x-1,y+1)];

            // Laplacian, force 1, spreading mass.
            let lap_u = 0.5*(u_rt + u_lt + u_dn + u_up) + 0.25*(u_ur + u_ul + u_dr + u_dl) - 3.0*u_ct;
            let lap_v = 0.5*(v_rt + v_lt + v_dn + v_up) + 0.25*(v_ur + v_ul + v_dr + v_dl) - 3.0*v_ct;
            let lap_w = 0.5*(w_rt + w_lt + w_dn + w_up) + 0.25*(w_ur + w_ul + w_dr + w_dl) - 3.0*w_ct;

            let mut force_u = lap_u;
            let mut force_v = lap_v;
            let mut force_w = lap_w;

            // Vacuum potential, force 2, pull w towards -1 aka vacuum state.
            let mass_squared = 0.02; 
            force_w -= mass_squared * (1.0 + w_ct);

            // 4. THE SKYRME PRESSURE (Force C)
            // Calculate the spatial gradients
            let du_dx = (u_rt - u_lt) * 0.5; let du_dy = (u_dn - u_up) * 0.5;
            let dv_dx = (v_rt - v_lt) * 0.5; let dv_dy = (v_dn - v_up) * 0.5;
            let dw_dx = (w_rt - w_lt) * 0.5; let dw_dy = (w_dn - w_up) * 0.5;

            let cross_u = dv_dx * dw_dy - dw_dx * dv_dy;
            let cross_v = dw_dx * du_dy - du_dx * dw_dy;
            let cross_w = du_dx * dv_dy - dv_dx * du_dy;

            let skyrme_coupling = 20.0; // up if skyrme shrinks, down if it explodes
            
            force_u += cross_u * skyrme_coupling;
            force_v += cross_v * skyrme_coupling;
            force_w += cross_w * skyrme_coupling;

            // advection, for local fluid momentum
            let mut vx = self.px[idx];
            let mut vy = self.py[idx];

            // speed of time limiter
            let max_v = 0.45; 
            let v_mag = (vx * vx + vy * vy).sqrt();
            if v_mag > max_v {
                vx = (vx / v_mag) * max_v;
                vy = (vy / v_mag) * max_v;
            }

            //advectioj force, then drag
            let advect_u = -(vx * du_dx + vy * du_dy);
            let advect_v = -(vx * dv_dx + vy * dv_dy);
            let advect_w = -(vx * dw_dx + vy * dw_dy);
            force_u += advect_u * 0.5;
            force_v += advect_v * 0.5;
            force_w += advect_w * 0.5;

            //momentum decay
            next_px[idx] = vx * 0.98;
            next_py[idx] = vy * 0.98;

            let speed = 0.12;

            let mut nu = (2.0 * u_ct) - u_pt + (speed * force_u);
            let mut nv = (2.0 * v_ct) - v_pt + (speed * force_v);
            let mut nw = (2.0 * w_ct) - w_pt + (speed * force_w);

            // 6. Hard Geometric Projection (Normalization)
            let mag = (nu*nu + nv*nv + nw*nw).sqrt();
            if mag > 1e-6 {
                nu /= mag; nv /= mag; nw /= mag;
            }

            // 7. Friction (Critical for long-term stability)
            // This "cools" the system so it settles into its ground state.
            let friction = 0.9995;
            nu = (nu - u_ct) * friction + u_ct;
            nv = (nv - v_ct) * friction + v_ct;
            nw = (nw - w_ct) * friction + w_ct;

            // 8. Final Save
            next_u[idx] = nu; next_v[idx] = nv; next_w[idx] = nw;
            next_u_prev[idx] = u_ct; next_v_prev[idx] = v_ct; next_w_prev[idx] = w_ct;
            
            // Matter visualization (Knot strength)
            let jacobian = (du_dx * dv_dy - du_dy * dv_dx).abs();
            next_r[idx] = (jacobian * 15.0).clamp(0.0, 1.0);
        }
    }
    self.u = next_u; self.u_prev = next_u_prev;
    self.v = next_v; self.v_prev = next_v_prev;
    self.w = next_w; self.w_prev = next_w_prev;
    self.rho = next_r;
    self.charge = self.recompute_charge_from_uv(&self.u, &self.v);
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
        (255 << 24) | (c << 16)
    } else if q < -0.02 {
        (255 << 24) | c
    } else {
        let g = (40.0 * rho.clamp(0.0, 1.0)) as u32;
        (255 << 24) | (g << 16) | (g << 8) | g
    }
}

fn core_to_color(core: f32) -> u32 {
    let c = (255.0 * core.clamp(0.0, 1.0)) as u32;
    (255 << 24) | (c << 16) | (c << 8) | c
}

fn skyrmion_color(u: f32, v: f32, w: f32) -> u32 {
    // Phase dictates the color (the U/V equator twist)
    let theta = v.atan2(u);
    
    // W dictates the brightness! 
    // Vacuum (W=-1) becomes 0.0 (Black). Core (W=1) becomes 1.0 (Bright).
    let brightness = (w * 0.5 + 0.5).clamp(0.0, 1.0);

    let red = ((theta.sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;
    let grn = (((theta + 2.0).sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;
    let blu = (((theta + 4.0).sin() * 0.5 + 0.5) * 255.0 * brightness) as u32;

    (255 << 24) | (red << 16) | (grn << 8) | blu
}

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer = vec![0u32; WINDOW_WIDTH * WINDOW_HEIGHT];
    let mut stats_frame: u32 = 0;
    let mut last_title = String::from("Starting...");

    let mut window = Window::new(
        "Now Wave | LMB + | RMB - | Shift+L push right | Shift+R push left | H heat | T warp | P rho | C charge | Space reset",
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
                        universe.charge = universe.recompute_charge_from_uv(&universe.u, &universe.v);
                    }
                    if is_r {
                        universe.spawn_vortex(cx, cy, -1.0);
                        universe.charge = universe.recompute_charge_from_uv(&universe.u, &universe.v);
                    }
                }
            }
        }

        last_m = is_l || is_r;

        let steps = if window.is_key_down(Key::T) { 100 } else { 5 };

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
                // We don't strictly need to recompute charge here, 
                // the step function handles it smoothly now!
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
        let d = ((peaks[0].0 - peaks[1].0).powi(2) + (peaks[0].1 - peaks[1].1).powi(2)).sqrt();
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

        let show_core = window.is_key_down(Key::X);
        let show_potential = window.is_key_down(Key::P);
        let show_charge = window.is_key_down(Key::C);
        let show_force = window.is_key_down(Key::F);

        for i in 0..buffer.len() {
            if show_core {
                let x = (i % SIM_WIDTH) as isize;
                let y = (i / SIM_WIDTH) as isize;
                let core = universe.local_core_strength(x, y);
                buffer[i] = core_to_color(core);
            } else if show_force { 
                // 2. Add the Diagnostic Heatmap
                // Calculate the absolute magnitude of the fluid momentum
                let mom = (universe.px[i].powi(2) + universe.py[i].powi(2)).sqrt();
                
                // Multiply by 20.0 to make even microscopic forces visibly glow red
                let heat = (mom * 20.0).clamp(0.0, 1.0);
                let c = (heat * 255.0) as u32;
                
                // Render as pure red on a black background
                buffer[i] = (255 << 24) | (c << 16) | (0 << 8) | 0;
            } else if show_charge {
                buffer[i] = charge_to_color(universe.charge[i], universe.rho[i]);
            } else {
                // Render the true 3D Skyrmion!
                buffer[i] = skyrmion_color(universe.u[i], universe.v[i], universe.w[i]);
            }
        }
        window
            .update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)
            .unwrap();
    }
}
