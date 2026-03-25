use minifb::{Key, Window, WindowOptions};

const SIM_WIDTH: usize = 200;
const SIM_HEIGHT: usize = 200;
const SCALE: usize = 4;
const WINDOW_WIDTH: usize = SIM_WIDTH * SCALE;
const WINDOW_HEIGHT: usize = SIM_HEIGHT * SCALE;

pub struct Universe {
    width: usize,
    height: usize,
    current: Vec<[f32; 2]>,
    previous: Vec<[f32; 2]>,
}

impl Universe {
    pub fn new(width: usize, height: usize) -> Self {
        let size = width * height;
        let current = vec![[0.0, 1.0]; size];  // Flat zero
        let previous = vec![[0.0, 1.0]; size]; // Flat zero

        Self { width, height, current, previous }
    }

    fn get_index(&self, x: isize, y: isize) -> usize {
        let wrapped_x = x.rem_euclid(self.width as isize) as usize;
        let wrapped_y = y.rem_euclid(self.height as isize) as usize;
        wrapped_y * self.width + wrapped_x
    }

pub fn step(&mut self) {
        let mut next = vec![[0.0, 0.0]; self.width * self.height];
        let wave_speed = 0.15; // Lower for vector stability

        for y in 0..self.height as isize {
            for x in 0..self.width as isize {
                let idx = self.get_index(x, y);

                // Get neighbor vectors
                let up = self.current[self.get_index(x, y - 1)];
                let down = self.current[self.get_index(x, y + 1)];
                let left = self.current[self.get_index(x - 1, y)];
                let right = self.current[self.get_index(x + 1, y)];

                let curr = self.current[idx];
                let prev = self.previous[idx];

                // Calculate Vector Laplacian (Difference for both U and V)
                let lap_u = up[0] + down[0] + left[0] + right[0] - (4.0 * curr[0]);
                let lap_v = up[1] + down[1] + left[1] + right[1] - (4.0 * curr[1]);


                // Update each component
                let mut next_u = (2.0 * curr[0]) - prev[0] + (wave_speed * lap_u);
                let mut next_v = (2.0 * curr[1]) - prev[1] + (wave_speed * lap_v);

                let potential_u = (0.0 - next_u) * 0.004;
                let potential_v = (1.0 - next_v) * 0.004;

                next_u += potential_u;
                next_v += potential_v;

                // RELATIVITY FIX: Normalize the vector
                // This ensures the "Time Needle" always has a length of 1.0. 
                // It can rotate, but it cannot grow or shrink.
                let mag = (next_u * next_u + next_v * next_v).sqrt();
                if mag > 0.0 {
                    next_u /= mag;
                    next_v /= mag;
                }

                next[idx] = [next_u, next_v];
            }
        }
        self.previous = self.current.clone();
        self.current = next;
    }
}

// Helper function to turn a float phase into an ARGB pixel
fn phase_to_color(vec: [f32; 2]) -> u32 {
    // Calculate the angle of the vector (atan2 gives -PI to PI)
    let angle = vec[1].atan2(vec[0]);
    
    // Map angle to RGB
    let r = ((angle.sin() * 0.5 + 0.5) * 255.0) as u32;
    let g = ((angle.cos() * 0.5 + 0.5) * 255.0) as u32;
    let b = (((angle + 2.0).sin() * 0.5 + 0.5) * 255.0) as u32;

    (255 << 24) | (r << 16) | (g << 8) | b
}

fn main() {
    let mut universe = Universe::new(SIM_WIDTH, SIM_HEIGHT);
    let mut buffer: Vec<u32> = vec![0; WINDOW_WIDTH * WINDOW_HEIGHT];

    let mut window = Window::new(
        "Now Wave Simulator - Sandbox Mode",
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        WindowOptions::default(),
    ).unwrap_or_else(|e| {
        panic!("{}", e);
    });

    window.set_target_fps(60);

    let mut was_mouse_down = false;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        
        let is_mouse_down = window.get_mouse_down(minifb::MouseButton::Left);

        // ONLY fire on the exact frame the mouse is clicked
         if is_mouse_down && !was_mouse_down {
            if let Some((mouse_x, mouse_y)) = window.get_mouse_pos(minifb::MouseMode::Discard) {
                let center_x = (mouse_x as usize / SCALE) as isize;
                let center_y = (mouse_y as usize / SCALE) as isize;
                let radius = 40.0_f32;
                let core_scale = 22.0_f32;

                for dy in -(radius as isize)..=(radius as isize) {
                    for dx in -(radius as isize)..=(radius as isize) {
                        let r2 = (dx * dx + dy * dy) as f32;
                        if r2 <= radius * radius {
                            let x_idx = (center_x + dx).rem_euclid(SIM_WIDTH as isize) as usize;
                            let y_idx = (center_y + dy).rem_euclid(SIM_HEIGHT as isize) as usize;
                            let idx = y_idx * SIM_WIDTH + x_idx;

                            let r = r2.sqrt();

                            // Tangential swirl direction
                            let mut tu = -(dy as f32);
                            let mut tv = dx as f32;
                            let tmag = (tu * tu + tv * tv).sqrt();

                            if tmag > 0.0 {
                                tu /= tmag;
                                tv /= tmag;
                            } else {
                                tu = 0.0;
                                tv = 1.0;
                            }

                            // Strength fades with radius
                            let twist = (-(r * r) / (2.0 * core_scale * core_scale)).exp();

                            // Blend between flat background [0,1] and vortex direction
                            let mut u = (1.0 - twist) * 0.0 + twist * tu;
                            let mut v = (1.0 - twist) * 1.0 + twist * tv;

                            let mag = (u * u + v * v).sqrt();
                            if mag > 0.0 {
                                u /= mag;
                                v /= mag;
                            }

                            universe.current[idx] = [u, v];
                            universe.previous[idx] = [u, v];
                        }
                    }
                }
            }
        }       
        // Update the mouse state for the next frame
        was_mouse_down = is_mouse_down;

        // 1. Advance the physics
        universe.step();

        // 2. Render colors
        for sy in 0..SIM_HEIGHT {
            for sx in 0..SIM_WIDTH {
                let color = phase_to_color(universe.current[sy * SIM_WIDTH + sx]);

                for oy in 0..SCALE {
                    for ox in 0..SCALE {
                        let wx = sx * SCALE + ox;
                        let wy = sy * SCALE + oy;
                        buffer[wy * WINDOW_WIDTH + wx] = color;
                    }
                }
            }
        }

        // 3. Draw
        window
            .update_with_buffer(&buffer, WINDOW_WIDTH, WINDOW_HEIGHT)
            .unwrap();
    }
}
