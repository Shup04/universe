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

fn spawn_vortex(&mut self, cx: isize, cy: isize, sign: f32) {
    let w = self.width as isize;
    let h = self.height as isize;
    
    let separation = 18.0; 
    let core_size = 4.0;

    // THE KICKSTART: We inject a vertical velocity.
    // A horizontal pair naturally propels itself UP (-Y) across the screen.
    let vy = -1.0_f32; 

    for y in 0..h {
        for x in 0..w {
            let idx = self.get_index(x, y);

            let mut dx = x - cx;
            let mut dy = y - cy;
            if dx > w/2 { dx -= w; } else if dx < -w/2 { dx += w; }
            if dy > h/2 { dy -= h; } else if dy < -h/2 { dy += h; }
            
            // --- 1. CURRENT FRAME ---
            let fx1 = dx as f32 + separation;
            let fy1 = dy as f32;
            let r1_sq = fx1 * fx1 + fy1 * fy1;

            let fx2 = dx as f32 - separation;
            let fy2 = dy as f32;
            let r2_sq = fx2 * fx2 + fy2 * fy2;

            let amp1 = (r1_sq.sqrt() / core_size).min(1.0);
            let amp2 = (r2_sq.sqrt() / core_size).min(1.0);
            let amp = amp1 * amp2;

            let theta1 = sign * fy1.atan2(fx1);
            let theta2 = -sign * fy2.atan2(fx2);
            let total_theta = theta1 + theta2 + std::f32::consts::PI / 2.0;

            let knot_u = amp * total_theta.cos();
            let knot_v = amp * total_theta.sin();

            // --- 2. PREVIOUS FRAME (The Time Shift) ---
            // If the particle is moving UP (-Y), then in the past, it was LOWER (+Y).
            // We shift the Y coordinates backward by 'vy' to create instant momentum.
            let fy1_prev = fy1 - vy; 
            let fy2_prev = fy2 - vy;

            let theta1_prev = sign * fy1_prev.atan2(fx1);
            let theta2_prev = -sign * fy2_prev.atan2(fx2);
            let total_theta_prev = theta1_prev + theta2_prev + std::f32::consts::PI / 2.0;

            let prev_knot_u = amp * total_theta_prev.cos();
            let prev_knot_v = amp * total_theta_prev.sin();

            // 3. Superposition onto the vacuum
            let old_u = self.u[idx];
            let old_v = self.v[idx];
            self.u[idx] = old_u * knot_u - old_v * knot_v;
            self.v[idx] = old_u * knot_v + old_v * knot_u;
            self.u_prev[idx] = self.u[idx];
            self.v_prev[idx] = self.v[idx];
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

            // --- 1. THE EMERGENCY BRAKES (Fixing the Coasting Bug) ---
            // Single particles max out at ~0.4 density. 
            // When they crash, stress spikes the density to 1.0. 
            // We only freeze when density hits 0.8, strictly catching collisions!
            let is_frozen = false;

            // --- 1. CALCULATE WAVE TENSION FIRST ---
            // We move this up so we can use it to detect the collision!
            let diff_up = (u_up - u_ct).powi(2) + (v_up - v_ct).powi(2);
            let diff_dn = (u_dn - u_ct).powi(2) + (v_dn - v_ct).powi(2);
            let diff_lt = (u_lt - u_ct).powi(2) + (v_lt - v_ct).powi(2);
            let diff_rt = (u_rt - u_ct).powi(2) + (v_rt - v_ct).powi(2);
            let phase_variance = diff_up + diff_dn + diff_lt + diff_rt;

            // --- 1. QUADRATIC TIME DILATION (The Honey Trap) ---
            // Squaring `cr` means time flows normally on the edges, 
            // but grinds to a 1% crawl perfectly at the dense core.
            let collision_spike = (phase_variance - 1.2).max(0.0);
            let local_speed = base_speed * (1.0 - 0.15 * cr).clamp(0.85, 1.0);

            // --- THE WAVE EQUATION ---
            let lap_u = u_up + u_dn + u_lt + u_rt - (4.0 * u_ct);
            let lap_v = v_up + v_dn + v_lt + v_rt - (4.0 * v_ct);

            // HIT THE BRAKES: By replacing the 'past' frame with the 'current' frame, 
            // we mathematically delete their momentum so they can't coast through the freeze.
            let eff_u_pt = if is_frozen { u_ct } else { u_pt };
            let eff_v_pt = if is_frozen { v_ct } else { v_pt };

            let mut nu = (2.0 * u_ct) - eff_u_pt + (local_speed * lap_u);
            let mut nv = (2.0 * v_ct) - eff_v_pt + (local_speed * lap_v);

            // --- THE MEXICAN HAT POTENTIAL ---
            let amp_sq = u_ct * u_ct + v_ct * v_ct;
            let restoring_force = 1.0 - amp_sq;
            nu += u_ct * restoring_force * 0.1;
            nv += v_ct * restoring_force * 0.1;

            // Universal friction
            nu = (nu - u_ct) * 0.997 + u_ct;
            nv = (nv - v_ct) * 0.997 + v_ct;

            // --- KINETIC ENERGY & TOPOLOGY ---
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

            // --- 2. VACUUM STRESS ---
            let stress = restoring_force.max(0.0);
            let target_r = (0.45 * phase_variance + 0.35 * knot_strength + 0.20 * stress).clamp(0.0, 1.0);

            // --- YOUR ADVECTION SOLVER (The Strong Force) ---
            let rho_up = self.rho[self.get_index(x, y-1)];
            let rho_dn = self.rho[self.get_index(x, y+1)];
            let rho_lt = self.rho[self.get_index(x-1, y)];
            let rho_rt = self.rho[self.get_index(x+1, y)];

            let strong_repulsion = 0.0; 
            let core_fx = (rho_lt.powi(3) - rho_rt.powi(3)) * strong_repulsion;
            let core_fy = (rho_up.powi(3) - rho_dn.powi(3)) * strong_repulsion;

            let fx = (rho_lt - rho_rt) * 0.003 + core_fx;
            let fy = (rho_up - rho_dn) * 0.003 + core_fy;

            //redo fx and fx to work with charge

            let q_up = self.charge[self.get_index(x, y - 1)];
            let q_dn = self.charge[self.get_index(x, y + 1)];
            let q_lt = self.charge[self.get_index(x - 1, y)];
            let q_rt = self.charge[self.get_index(x + 1, y)];

            let q0 = self.charge[idx];
            let bind = 0.006;

            let bind_fx = if q0 > 0.0 {
                (q_rt.min(0.0).abs() - q_lt.min(0.0).abs()) * bind
            } else {
                (q_lt.max(0.0) - q_rt.max(0.0)) * bind
            };

            let bind_fy = if q0 > 0.0 {
                (q_dn.min(0.0).abs() - q_up.min(0.0).abs()) * bind
            } else {
                (q_up.max(0.0) - q_dn.max(0.0)) * bind
            };

            let fx = fx + bind_fx;
            let fy = fy + bind_fy;

                //-----

            let pxi = self.px[idx];
            let pyi = self.py[idx];
            let drag = 0.995;
            let mut npx = pxi * drag + fx;
            let mut npy = pyi * drag + fy;

            let src_x = x as f32 - npx;
            let src_y = y as f32 - npy;
            let src_idx = self.get_index(src_x.floor() as isize, src_y.floor() as isize);

            // --- 3. CLEAN MATTER DECAY (No more Ice-Nine crystals) ---
            let advected_rho = self.rho[src_idx];
            let formed_rho = (cr * 0.92) + (target_r * 0.08);
            let mut nr = (0.75 * formed_rho + 0.25 * advected_rho).clamp(0.001, 1.0);

            // --- 2. THE TOTAL FREEZE ---
            //if is_frozen {
             //   nr = cr.max(0.85); // Lock the matter so it never evaporates back down
             //   npx = 0.0;         // Kill the fluid pressure
            //    npy = 0.0;
            //}

            // --- SAVE TO NEXT FRAME ---
            next_u[idx] = nu;          
            next_u_prev[idx] = u_ct;   
            next_v[idx] = nv;          
            next_v_prev[idx] = v_ct;   
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

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer = vec![0u32; WINDOW_WIDTH * WINDOW_HEIGHT];

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

        if peaks.len() == 2 {
            let d = ((peaks[0].0 - peaks[1].0).powi(2) + (peaks[0].1 - peaks[1].1).powi(2)).sqrt();
            window.set_title(&format!(
                "Dist: {:.4} | Total charge: {:.4} | Total momentum: {:.4}",
                d, total_charge, total_momentum
            ));
        } else {
            window.set_title(&format!(
                "Total charge: {:.4} | Total momentum: {:.4}",
                total_charge, total_momentum
            ));
        }

        let show_potential = window.is_key_down(Key::P);
        let show_charge = window.is_key_down(Key::C);

        for i in 0..buffer.len() {
            // Extract the angle on-the-fly for the renderer
            let current_theta = universe.v[i].atan2(universe.u[i]);
            
            buffer[i] = if show_charge {
                charge_to_color(universe.charge[i], universe.rho[i])
            } else {
                theta_rho_to_color(current_theta, universe.rho[i], show_potential)
            };
        }

        window
            .update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)
            .unwrap();
    }
}
