use minifb::{Key, Window, WindowOptions};
use std::f32::consts::PI;
use rand::Rng;

const SIM_WIDTH: usize = 400; // Scaled down for speed, set to 800 if preferred
const SIM_HEIGHT: usize = 400;
const SCALE: usize = 2; 
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
        
        let base_speed = 0.22; 
        let latch = 0.45;            
        let g_strength = 15.0;       
        let surface_tension = 0.2; 
        let friction = 0.95;        

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

                let local_speed = base_speed / (1.0 + (cr * cr) * g_strength).min(10.0);
                let lap = (up-ct).sin() + (dn-ct).sin() + (lt-ct).sin() + (rt-ct).sin();
                let mut nt = (2.0 * ct) - pt + (local_speed * lap);

                // --- THE GRAVITY HOOK ---
                // This makes the phase "fall" into areas of high density.
                // It creates a 'self-reinforcing' bridge.
                nt += (cr * 0.05);
                nt -= ct.sin() * latch;
                
                if cr < 0.1 { nt = (nt - ct) * friction + ct; }

                let grad = (up-dn).powi(2) + (lt-rt).powi(2);
                let target_r = if grad > surface_tension { (grad * 3.0).min(1.0) } else { 0.0 };
                let nr = (cr * 0.95) + (target_r * 0.05);

                next_t[idx] = nt;
                next_r[idx] = nr;
            }
        }
        self.theta_prev = self.theta.clone();
        self.theta = next_t;
        self.rho = next_r;
    }
}

fn phase_to_color(theta: f32, rho: f32, show_potential: bool) -> u32 {
    let r = (theta.sin() * 0.5 + 0.5) * 255.0;
    let g = ((theta + 2.0).sin() * 0.5 + 0.5) * 255.0;
    let b = ((theta + 4.0).sin() * 0.5 + 0.5) * 255.0;

    let brightness = if show_potential {
        (rho * 4.0).clamp(0.05, 1.0) // Potential view makes the vacuum "glow"
    } else {
        if rho < 0.05 { 0.0 } else { rho.clamp(0.0, 1.0) }
    };
    
    let red = (r * brightness) as u32;
    let green = (g * brightness) as u32;
    let blue = (b * brightness) as u32;

    (255 << 24) | (red << 16) | (green << 8) | blue
}

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer = vec![0u32; WINDOW_WIDTH * WINDOW_HEIGHT];
    let mut window = Window::new("Now Wave - [H] Heat [P] Potential [N] Nudge", WINDOW_WIDTH, WINDOW_HEIGHT, WindowOptions::default()).unwrap();
    let mut last_m = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        if window.is_key_down(Key::Space) { universe = Universe::new(SIM_WIDTH, SIM_HEIGHT); }
        
        // 1. INPUT HANDLING
        let is_m = window.get_mouse_down(minifb::MouseButton::Left) || window.get_mouse_down(minifb::MouseButton::Right);
        if is_m && !last_m {
            if let Some((mx, my)) = window.get_mouse_pos(minifb::MouseMode::Discard) {
                let (cx, cy) = (mx as isize / SCALE as isize, my as isize / SCALE as isize);
                let swirl = if window.get_mouse_down(minifb::MouseButton::Right) { -1.0 } else { 1.0 };
                for dy in -20..20 {
                    for dx in -20..20 {
                        if dx*dx + dy*dy <= 400 {
                            let idx = universe.get_index(cx+dx, cy+dy);
                            let angle = (dy as f32).atan2(dx as f32);
                            let intensity = (-( (dx*dx+dy*dy) as f32 ) / 100.0).exp();
                            universe.theta[idx] = angle * swirl * intensity;
                            universe.theta_prev[idx] = universe.theta[idx];
                            universe.rho[idx] = 0.8 * intensity;
                        }
                    }
                }
            }
        }
        last_m = is_m;

        // 2. INTERACTIVE TOOLS
        if window.is_key_down(Key::H) { // HEAT: Jitter all matter
            let mut rng = rand::thread_rng();
            for i in 0..universe.theta.len() {
                if universe.rho[i] > 0.1 { universe.theta[i] += rng.gen_range(-0.1..0.1); }
            }
        }
        if window.is_key_down(Key::N) { // NUDGE: Small momentum push
            let mut rng = rand::thread_rng();
            for i in 0..universe.theta.len() {
                if universe.rho[i] > 0.1 { universe.theta_prev[i] += rng.gen_range(-0.02..0.02); }
            }
        }
        let show_potential = window.is_key_down(Key::P);

        universe.step();

        // 3. HUD LOGIC
        // --- SUB-PIXEL CENTROID HUD ---
        let mut centers = Vec::new();
        let mut search_rho = universe.rho.clone();

        for _ in 0..2 {
            let mut total_rho = 0.0;
            let mut sum_x = 0.0;
            let mut sum_y = 0.0;
            let mut max_val = 0.0;
            let mut max_idx = 0;

            // Find the core of the next particle
            for (i, &val) in search_rho.iter().enumerate() {
                if val > max_val { max_val = val; max_idx = i; }
            }

            if max_val > 0.3 {
                let px = (max_idx % SIM_WIDTH) as isize;
                let py = (max_idx / SIM_WIDTH) as isize;
                
                // Calculate the "Center of Mass" in a 20px radius
                for dy in -20..20 {
                    for dx in -20..20 {
                        let idx = universe.get_index(px + dx, py + dy);
                        let val = search_rho[idx];
                        if val > 0.1 {
                            sum_x += (px + dx) as f32 * val;
                            sum_y += (py + dy) as f32 * val;
                            total_rho += val;
                            search_rho[idx] = 0.0; // Clear for next particle
                        }
                    }
                }
                centers.push((sum_x / total_rho, sum_y / total_rho));
            }
        }

        if centers.len() == 2 {
            let (c1, c2) = (centers[0], centers[1]);
            let dx = c2.0 - c1.0;
            let dy = c2.1 - c1.1;
            let dist = (dx*dx + dy*dy).sqrt();
            window.set_title(&format!("Sub-Pixel Dist: {:.4} | Force: {:.5}", dist, 0.6)); // 0.6 is your stable bridge
        }

        // 4. RENDER
        for y in 0..SIM_HEIGHT {
            for x in 0..SIM_WIDTH {
                let idx = y * SIM_WIDTH + x;
                let color = phase_to_color(universe.theta[idx], universe.rho[idx], show_potential);
                for oy in 0..SCALE {
                    for ox in 0..SCALE {
                        buffer[(y * SCALE + oy) * WINDOW_WIDTH + (x * SCALE + ox)] = color;
                    }
                }
            }
        }
        window.update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT).unwrap();
    }
}
