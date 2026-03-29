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
            u: vec![1.0; size],
            u_prev: vec![1.0; size],
            v: vec![0.0; size],
            v_prev: vec![0.0; size],
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

                let c_theta  = v_field[self.get_index(x, y)].atan2(u_field[self.get_index(x, y)]);
                let r_theta  = v_field[self.get_index(x + 1, y)].atan2(u_field[self.get_index(x + 1, y)]);
                let ur_theta = v_field[self.get_index(x + 1, y - 1)].atan2(u_field[self.get_index(x + 1, y - 1)]);
                let u_theta  = v_field[self.get_index(x, y - 1)].atan2(u_field[self.get_index(x, y - 1)]);

                let d1 = wrap_angle(r_theta - c_theta);
                let d2 = wrap_angle(ur_theta - r_theta);
                let d3 = wrap_angle(u_theta - ur_theta);
                let d4 = wrap_angle(c_theta - u_theta);

                next_charge[idx] = (d1 + d2 + d3 + d4) / (2.0 * PI);
            }
        }

        next_charge
    }

    fn amp_at(&self, idx: usize) -> f32 {
        (self.u[idx] * self.u[idx] + self.v[idx] * self.v[idx]).sqrt()
    }

    fn local_phase_variance(&self, x: isize, y: isize) -> f32 {
        let idx = self.get_index(x, y);

        let u_ct = self.u[idx];
        let v_ct = self.v[idx];

        let u_up = self.u[self.get_index(x, y - 1)];
        let u_dn = self.u[self.get_index(x, y + 1)];
        let u_lt = self.u[self.get_index(x - 1, y)];
        let u_rt = self.u[self.get_index(x + 1, y)];

        let v_up = self.v[self.get_index(x, y - 1)];
        let v_dn = self.v[self.get_index(x, y + 1)];
        let v_lt = self.v[self.get_index(x - 1, y)];
        let v_rt = self.v[self.get_index(x + 1, y)];

        let diff_up = (u_up - u_ct).powi(2) + (v_up - v_ct).powi(2);
        let diff_dn = (u_dn - u_ct).powi(2) + (v_dn - v_ct).powi(2);
        let diff_lt = (u_lt - u_ct).powi(2) + (v_lt - v_ct).powi(2);
        let diff_rt = (u_rt - u_ct).powi(2) + (v_rt - v_ct).powi(2);

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
    
    // The Soliton Drive configuration
    let separation = 10.0; 
    let core_size = 4.0;

    // THE KICKSTART: This gives them the momentum to create 
    // the hydrodynamic pressure wake that stops them from spinning!
    let vy = -0.5_f32; 

    for y in 0..h {
        for x in 0..w {
            let idx = self.get_index(x, y);

            let mut dx = x - cx;
            let mut dy = y - cy;
            if dx > w/2 { dx -= w; } else if dx < -w/2 { dx += w; }
            if dy > h/2 { dy -= h; } else if dy < -h/2 { dy += h; }
            
            // Core 1 (+1)
            let fx1 = dx as f32 + separation;
            let fy1 = dy as f32;
            let r1_sq = fx1 * fx1 + fy1 * fy1;

            // Core 2 (-1)
            let fx2 = dx as f32 - separation;
            let fy2 = dy as f32;
            let r2_sq = fx2 * fx2 + fy2 * fy2;

            let amp1 = (r1_sq.sqrt() / core_size).min(1.0);
            let amp2 = (r2_sq.sqrt() / core_size).min(1.0);
            let amp = amp1 * amp2;

            let theta1 = 1.0 * fy1.atan2(fx1);
            let theta2 = -1.0 * fy2.atan2(fx2);
            let total_theta = theta1 + theta2 + std::f32::consts::PI / 2.0;

            let knot_u = amp * total_theta.cos();
            let knot_v = amp * total_theta.sin();

            // The Time Shift (Applying the forward momentum)
            let fy1_prev = fy1 - vy; 
            let fy2_prev = fy2 - vy;
            let theta1_prev = 1.0 * fy1_prev.atan2(fx1);
            let theta2_prev = -1.0 * fy2_prev.atan2(fx2);
            let total_theta_prev = theta1_prev + theta2_prev + std::f32::consts::PI / 2.0;

            let prev_knot_u = amp * total_theta_prev.cos();
            let prev_knot_v = amp * total_theta_prev.sin();

            // Superposition onto the vacuum
            let old_u = self.u[idx];
            let old_v = self.v[idx];
            
            self.u[idx] = old_u * knot_u - old_v * knot_v;
            self.v[idx] = old_u * knot_v + old_v * knot_u;
            
            // Anchor the past frame shifted to prevent the rubber band snap!
            let old_u_prev = self.u_prev[idx];
            let old_v_prev = self.v_prev[idx];
            self.u_prev[idx] = old_u_prev * prev_knot_u - old_v_prev * prev_knot_v;
            self.v_prev[idx] = old_u_prev * prev_knot_v + old_v_prev * prev_knot_u;
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
    // 1. Create the empty vectors for the NEXT frame
    let mut next_u = vec![0.0; self.u.len()];
    let mut next_u_prev = vec![0.0; self.u.len()];
    let mut next_v = vec![0.0; self.v.len()];
    let mut next_v_prev = vec![0.0; self.v.len()];
    let mut next_r = vec![0.0; self.rho.len()];
    let mut next_px = vec![0.0; self.px.len()];
    let mut next_py = vec![0.0; self.py.len()];
    
    let base_speed = 0.4;

    for y in 0..self.height as isize {
        for x in 0..self.width as isize {
            let idx = self.get_index(x, y);

            // --- GATHER SPATIAL & TIME DATA FOR U ---
            let u_up = self.u[self.get_index(x, y-1)];
            let u_dn = self.u[self.get_index(x, y+1)];
            let u_lt = self.u[self.get_index(x-1, y)];
            let u_rt = self.u[self.get_index(x+1, y)];
            let u_ct = self.u[idx];
            let u_pt = self.u_prev[idx];

            // --- GATHER SPATIAL & TIME DATA FOR V ---
            let v_up = self.v[self.get_index(x, y-1)];
            let v_dn = self.v[self.get_index(x, y+1)];
            let v_lt = self.v[self.get_index(x-1, y)];
            let v_rt = self.v[self.get_index(x+1, y)];
            let v_ct = self.v[idx];
            let v_pt = self.v_prev[idx];

            let cr = self.rho[idx];

            // --- 1. TIME DILATION ---
            let local_speed = base_speed * (1.0 - 0.15 * cr).clamp(0.85, 1.0);

            // --- 2. THE WAVE EQUATION ---
            let lap_u = u_up + u_dn + u_lt + u_rt - (4.0 * u_ct);
            let lap_v = v_up + v_dn + v_lt + v_rt - (4.0 * v_ct);

            let mut nu = (2.0 * u_ct) - u_pt + (local_speed * lap_u);
            let mut nv = (2.0 * v_ct) - v_pt + (local_speed * lap_v);

            // --- 3. THE MEXICAN HAT POTENTIAL ---
            let amp_sq = u_ct * u_ct + v_ct * v_ct;
            let restoring_force = 1.0 - amp_sq;
            nu += u_ct * restoring_force * 0.1;
            nv += v_ct * restoring_force * 0.1;

            // Universal friction
            nu = (nu - u_ct) * 0.997 + u_ct;
            nv = (nv - v_ct) * 0.997 + v_ct;

            // --- 4. KINETIC ENERGY & TOPOLOGY ---
            let diff_up = (u_up - u_ct).powi(2) + (v_up - v_ct).powi(2);
            let diff_dn = (u_dn - u_ct).powi(2) + (v_dn - v_ct).powi(2);
            let diff_lt = (u_lt - u_ct).powi(2) + (v_lt - v_ct).powi(2);
            let diff_rt = (u_rt - u_ct).powi(2) + (v_rt - v_ct).powi(2);
            let phase_variance = diff_up + diff_dn + diff_lt + diff_rt;

            let du_dx = (u_rt - u_lt) * 0.5;
            let dv_dy = (v_dn - v_up) * 0.5; 
            let du_dy = (u_dn - u_up) * 0.5;
            let dv_dx = (v_rt - v_lt) * 0.5;

            let jacobian = du_dx * dv_dy - du_dy * dv_dx; 
            let knot_strength = (jacobian * 5.0).abs().min(1.0);

            let stress = restoring_force.max(0.0);
            let target_r = (0.45 * phase_variance + 0.35 * knot_strength + 0.20 * stress).clamp(0.0, 1.0);

            // --- 5. THE MACROSCOPIC RADAR (Fixing the Hurricane Eye) ---
            // By looking 3 pixels away, the dead-center of the core can "see" 
            // the massive collision front of the other particle and react to it!
            let rho_up3 = self.rho[self.get_index(x, y-3)];
            let rho_dn3 = self.rho[self.get_index(x, y+3)];
            let rho_lt3 = self.rho[self.get_index(x-3, y)];
            let rho_rt3 = self.rho[self.get_index(x+3, y)];

            // Because rho is naturally damped, this is unconditionally stable. No explosions.
            let cloud_repulsion = 2.0; 
            let fx = (rho_lt3.powi(3) - rho_rt3.powi(3)) * cloud_repulsion;
            let fy = (rho_up3.powi(3) - rho_dn3.powi(3)) * cloud_repulsion;

            let pxi = self.px[idx];
            let pyi = self.py[idx];
            let drag = 0.99; // Smooth, stable drag
            
            // --- 6. THE GOVERNOR ---
            let max_speed = 0.5; 
            let mut npx = (pxi * drag + fx).clamp(-max_speed, max_speed);
            let mut npy = (pyi * drag + fy).clamp(-max_speed, max_speed);

            // --- 7. THE DRIVE SHAFT ---
            let advect_u = (npx * du_dx + npy * du_dy).clamp(-0.20, 0.20);
            let advect_v = (npx * dv_dx + npy * dv_dy).clamp(-0.20, 0.20);

            nu -= advect_u;
            nv -= advect_v;
            
            let out_u_prev = u_ct - advect_u;
            let out_v_prev = v_ct - advect_v;

            let src_x = x as f32 - npx;
            let src_y = y as f32 - npy;
            let src_idx = self.get_index(src_x.floor() as isize, src_y.floor() as isize);

            // --- 8. CLEAN MATTER DECAY ---
            let advected_rho = self.rho[src_idx];
            let formed_rho = (cr * 0.94) + (target_r * 0.04);
            let nr = (0.82 * formed_rho + 0.18 * advected_rho).clamp(0.001, 1.0);

            // --- 9. SAVE TO NEXT FRAME ---
            next_u[idx] = nu;          
            next_u_prev[idx] = out_u_prev;
            next_v[idx] = nv;          
            next_v_prev[idx] = out_v_prev; 
            next_r[idx] = nr;
            next_px[idx] = npx;
            next_py[idx] = npy;
        } // End of x loop
    } // End of y loop
    
    // Lock it all in
    let new_charge = self.recompute_charge_from_uv(&next_u, &next_v);
    self.u_prev = next_u_prev;
    self.u = next_u;
    self.v_prev = next_v_prev;
    self.v = next_v;
    self.rho = next_r;
    self.px = next_px;
    self.py = next_py;
    self.charge = new_charge;
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
                    }big news recently
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
                let current_theta = universe.v[i].atan2(universe.u[i]);
                buffer[i] = theta_rho_to_color(current_theta, universe.rho[i], show_potential);
            }
        }
        window
            .update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)
            .unwrap();
    }
}
