use minifb::{Key, Window, WindowOptions, MouseMode};
use std::f32::consts::PI;
use rand::Rng;

const SIM_WIDTH: usize = 600;  // 800 might be slow in Debug mode, 600 is a good middle ground
const SIM_HEIGHT: usize = 600;
const SCALE: usize = 1; 
const WINDOW_WIDTH: usize = SIM_WIDTH * SCALE;
const WINDOW_HEIGHT: usize = SIM_HEIGHT * SCALE;

pub struct Universe {
    width: usize, height: usize,
    theta: Vec<f32>, theta_prev: Vec<f32>, rho: Vec<f32>,
}

impl Universe {
    pub fn new(width: usize, height: usize) -> Self {
        let size = width * height;
        Self { width, height, theta: vec![0.0; size], theta_prev: vec![0.0; size], rho: vec![0.01; size] }
    }

    fn get_index(&self, x: isize, y: isize) -> usize {
        let wx = x.rem_euclid(self.width as isize) as usize;
        let wy = y.rem_euclid(self.height as isize) as usize;
        wy * self.width + wx
    }

pub fn step(&mut self) {
    let mut next_t = vec![0.0; self.theta.len()];
    let mut next_r = vec![0.0; self.rho.len()];
    
    let base_speed = 0.4; 

    for y in 0..self.height as isize {
        for x in 0..self.width as isize {
            let idx = self.get_index(x, y);
            let up = self.theta[self.get_index(x, y-1)];
            let dn = self.theta[self.get_index(x, y+1)];
            let lt = self.theta[self.get_index(x-1, y)];
            let rt = self.theta[self.get_index(x+1, y)];

            let ct = self.theta[idx];
            let pt = self.theta_prev[idx];
            let cr = self.rho[idx];

            // 1. Gentle Time Dilation (NO FREEZING)
            // Time slows down by a maximum of 3x. This prevents the boundary explosion.
            let local_speed = base_speed / (1.0 + cr * 2.0);

            // 2. The Wave Equation
            let lap = (up-ct).sin() + (dn-ct).sin() + (lt-ct).sin() + (rt-ct).sin();
            let mut nt = (2.0 * ct) - pt + (local_speed * lap);

            // Universal friction to stop noise from echoing forever
            nt = (nt - ct) * 0.99 + ct;

            // --- 3. THE "PHASE VARIANCE" KNOT DETECTOR ---
            // We measure how violently the phase clashes with its neighbors.
            // At the exact core of a knot, this is high. At the edges, it's near zero.
            let diff_up = (up-ct).sin().powi(2);
            let diff_dn = (dn-ct).sin().powi(2);
            let diff_lt = (lt-ct).sin().powi(2);
            let diff_rt = (rt-ct).sin().powi(2);
            let phase_variance = diff_up + diff_dn + diff_lt + diff_rt;

            // 4. Density tied directly to the knot's core
            // Because variance is only high at the exact center, it CANNOT grow endlessly.
            let target_r = (phase_variance * 0.5).clamp(0.0, 1.0);
            let nr = (cr * 0.95) + (target_r * 0.05);

            next_t[idx] = nt;
            next_r[idx] = nr.max(0.01); // Safe vacuum floor
        }
    }
    self.theta_prev = self.theta.clone();
    self.theta = next_t;
    self.rho = next_r;
}
}

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer = vec![0u32; WINDOW_WIDTH * WINDOW_HEIGHT];
    let mut window = Window::new("Now Wave - T:Warp | H:Heat | P:Potential", WINDOW_WIDTH, WINDOW_HEIGHT, WindowOptions::default()).unwrap();
    let mut last_m = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if window.is_key_down(Key::Space) { universe = Universe::new(SIM_WIDTH, SIM_HEIGHT); }
        
        // --- 1. MOUSE INJECTION (DPI FIX) ---
let is_m = window.get_mouse_down(minifb::MouseButton::Left);
if is_m && !last_m {
    if let Some((mx, my)) = window.get_mouse_pos(minifb::MouseMode::Discard) {
        let (cx, cy) = (mx as isize / SCALE as isize, my as isize / SCALE as isize);
        
        // Inject a Dipole: Two vortices spinning opposite directions
        for dy in -40..40 {
            for dx in -60..60 {
                let idx = universe.get_index(cx+dx, cy+dy);
                
                // Left Knot (+1 Winding)
                let d1_sq = ((dx+20).pow(2) + dy.pow(2)) as f32;
                if d1_sq < 400.0 {
                    universe.theta[idx] += (dy as f32).atan2((dx+20) as f32);
                    universe.rho[idx] = 1.0;
                }
                
                // Right Knot (-1 Winding)
                let d2_sq = ((dx-20).pow(2) + dy.pow(2)) as f32;
                if d2_sq < 400.0 {
                    universe.theta[idx] -= (dy as f32).atan2((dx-20) as f32); // Notice the minus sign
                    universe.rho[idx] = 1.0;
                }
            }
        }
        universe.theta_prev = universe.theta.clone();
    }
}
        last_m = is_m;

        // --- 2. TIME WARP & HEAT ---
        let steps = if window.is_key_down(Key::T) { 100 } else { 1 };
        for _ in 0..steps {
            universe.step();
            if window.is_key_down(Key::H) {
                let mut rng = rand::thread_rng();
                for i in 0..universe.theta.len() {
                    if universe.rho[i] > 0.1 { universe.theta[i] += rng.gen_range(-0.1..0.1); }
                }
            }
        }

        // --- 3. CENTROID HUD (Sub-Pixel Tracking) ---
        let mut peaks = Vec::new();
        let mut search_rho = universe.rho.clone();
        for _ in 0..2 {
            let mut max_val = 0.0; let mut max_idx = 0;
            for (i, &val) in search_rho.iter().enumerate() { if val > max_val { max_val = val; max_idx = i; } }
            if max_val > 0.4 {
                let px = (max_idx % SIM_WIDTH) as isize;
                let py = (max_idx / SIM_WIDTH) as isize;
                let (mut sx, mut sy, mut total) = (0.0, 0.0, 0.0);
                for dy in -30..30 {
                    for dx in -30..30 {
                        let idx = universe.get_index(px+dx, py+dy);
                        let v = search_rho[idx];
                        if v > 0.1 { sx += (px+dx) as f32 * v; sy += (py+dy) as f32 * v; total += v; search_rho[idx] = 0.0; }
                    }
                }
                peaks.push((sx/total, sy/total));
            }
        }

        if peaks.len() == 2 {
            let d = ((peaks[0].0-peaks[1].0).powi(2) + (peaks[0].1-peaks[1].1).powi(2)).sqrt();
            window.set_title(&format!("Dist: {:.4}", d));
        }

        // --- 4. RENDER ---
        let p_view = window.is_key_down(Key::P);
        for i in 0..buffer.len() {
            let t = universe.theta[i];
            let r = universe.rho[i];
            let bri = if p_view { (r * 5.0).clamp(0.05, 1.0) } else { if r < 0.05 { 0.0 } else { r } };
            let red = ((t.sin() * 0.5 + 0.5) * 255.0 * bri) as u32;
            let grn = (((t+2.0).sin() * 0.5 + 0.5) * 255.0 * bri) as u32;
            let blu = (((t+4.0).sin() * 0.5 + 0.5) * 255.0 * bri) as u32;
            buffer[i] = (255 << 24) | (red << 16) | (grn << 8) | blu;
        }
        window.update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT).unwrap();
    }
}
