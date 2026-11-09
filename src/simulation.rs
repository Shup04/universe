use std::f32::consts::PI;

#[derive(Clone, Debug, PartialEq)]
pub struct Universe {
    pub(crate) width: usize,
    pub(crate) height: usize,
    pub(crate) rho: Vec<f32>,
    pub(crate) charge: Vec<f32>,
    pub(crate) px: Vec<f32>,
    pub(crate) py: Vec<f32>,
    pub(crate) u: Vec<f32>,
    pub(crate) u_prev: Vec<f32>,
    pub(crate) v: Vec<f32>,
    pub(crate) v_prev: Vec<f32>,
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
    pub fn reset(&mut self) {
        *self = Self::new(self.width, self.height);
    }

    pub fn is_finite(&self) -> bool {
        [
            &self.u,
            &self.v,
            &self.u_prev,
            &self.v_prev,
            &self.rho,
            &self.charge,
            &self.px,
            &self.py,
        ]
        .iter()
        .all(|field| field.iter().all(|value| value.is_finite()))
    }

    pub fn new(width: usize, height: usize) -> Self {
        assert!(width > 0 && height > 0, "grid dimensions must be positive");
        let size = width.checked_mul(height).expect("grid size overflow");
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

    pub(crate) fn get_index(&self, x: isize, y: isize) -> usize {
        let wx = x.rem_euclid(self.width as isize) as usize;
        let wy = y.rem_euclid(self.height as isize) as usize;
        wy * self.width + wx
    }

    pub(crate) fn recompute_charge_from_uv(&self, u_field: &[f32], v_field: &[f32]) -> Vec<f32> {
        let mut next_charge = vec![0.0; self.width * self.height];

        for y in 0..self.height as isize {
            for x in 0..self.width as isize {
                let idx = self.get_index(x, y);

                let c_theta = v_field[self.get_index(x, y)].atan2(u_field[self.get_index(x, y)]);
                let r_theta =
                    v_field[self.get_index(x + 1, y)].atan2(u_field[self.get_index(x + 1, y)]);
                let ur_theta = v_field[self.get_index(x + 1, y - 1)]
                    .atan2(u_field[self.get_index(x + 1, y - 1)]);
                let u_theta =
                    v_field[self.get_index(x, y - 1)].atan2(u_field[self.get_index(x, y - 1)]);

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

    pub(crate) fn local_core_strength(&self, x: isize, y: isize) -> f32 {
        let idx = self.get_index(x, y);
        let amp = self.amp_at(idx);
        let knot = self.local_knot_strength(x, y);
        let charge = self.charge[idx].abs().min(1.0);

        // core = strong topology + strong charge + low amplitude center
        (0.45 * knot + 0.45 * charge + 0.10 * (1.0 - amp)).clamp(0.0, 1.0)
    }

    #[cfg(feature = "gui")]
    pub(crate) fn local_pair_overlap(&self, x: isize, y: isize) -> f32 {
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

    /// Superpose a pair with opposite winding. Negative sign reverses both cores.
    /// Screen Y increases downward; the charge diagnostic walks upward around cells.
    pub fn spawn_vortex(&mut self, cx: isize, cy: isize, sign: f32) {
        assert!(sign == 1.0 || sign == -1.0, "orientation must be +1 or -1");
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
                if dx > w / 2 {
                    dx -= w;
                } else if dx < -w / 2 {
                    dx += w;
                }
                if dy > h / 2 {
                    dy -= h;
                } else if dy < -h / 2 {
                    dy += h;
                }

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

                let theta1 = sign * fy1.atan2(fx1);
                let theta2 = -sign * fy2.atan2(fx2);
                let total_theta = theta1 + theta2 + std::f32::consts::PI / 2.0;

                let knot_u = amp * total_theta.cos();
                let knot_v = amp * total_theta.sin();

                // The Time Shift (Applying the forward momentum)
                let fy1_prev = fy1 - vy;
                let fy2_prev = fy2 - vy;
                let theta1_prev = sign * fy1_prev.atan2(fx1);
                let theta2_prev = -sign * fy2_prev.atan2(fx2);
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
        self.charge = self.recompute_charge_from_uv(&self.u, &self.v);
    }

    #[cfg(any(feature = "gui", test))]
    pub(crate) fn push_nearby(&mut self, cx: isize, cy: isize, push_x: f32, push_y: f32) {
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
                let u_up = self.u[self.get_index(x, y - 1)];
                let u_dn = self.u[self.get_index(x, y + 1)];
                let u_lt = self.u[self.get_index(x - 1, y)];
                let u_rt = self.u[self.get_index(x + 1, y)];
                let u_ct = self.u[idx];
                let u_pt = self.u_prev[idx];

                // --- GATHER SPATIAL & TIME DATA FOR V ---
                let v_up = self.v[self.get_index(x, y - 1)];
                let v_dn = self.v[self.get_index(x, y + 1)];
                let v_lt = self.v[self.get_index(x - 1, y)];
                let v_rt = self.v[self.get_index(x + 1, y)];
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
                let target_r =
                    (0.45 * phase_variance + 0.35 * knot_strength + 0.20 * stress).clamp(0.0, 1.0);

                // --- 5. THE MACROSCOPIC RADAR (Fixing the Hurricane Eye) ---
                // By looking 3 pixels away, the dead-center of the core can "see"
                // the massive collision front of the other particle and react to it!
                let rho_up3 = self.rho[self.get_index(x, y - 3)];
                let rho_dn3 = self.rho[self.get_index(x, y + 3)];
                let rho_lt3 = self.rho[self.get_index(x - 3, y)];
                let rho_rt3 = self.rho[self.get_index(x + 3, y)];

                // Heuristic density-gradient repulsion; no unconditional stability claim.
                let cloud_repulsion = 2.0;
                let fx = (rho_lt3.powi(3) - rho_rt3.powi(3)) * cloud_repulsion;
                let fy = (rho_up3.powi(3) - rho_dn3.powi(3)) * cloud_repulsion;

                let pxi = self.px[idx];
                let pyi = self.py[idx];
                let drag = 0.99; // Smooth, stable drag

                // --- 6. THE GOVERNOR ---
                let max_speed = 0.5;
                let npx = (pxi * drag + fx).clamp(-max_speed, max_speed);
                let npy = (pyi * drag + fy).clamp(-max_speed, max_speed);

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

#[derive(Clone, Copy, Debug)]
pub enum View {
    Matter,
    Potential,
    Charge,
    Core,
    Momentum,
}

impl Universe {
    /// Packed RGB pixels shared by the window and headless captures.
    pub fn frame(&self, view: View) -> Vec<u32> {
        (0..self.u.len())
            .map(|i| match view {
                View::Core => core_to_color(
                    self.local_core_strength((i % self.width) as isize, (i / self.width) as isize),
                ),
                View::Momentum => {
                    let mom = (self.px[i].powi(2) + self.py[i].powi(2)).sqrt();
                    let c = ((mom * 20.0).clamp(0.0, 1.0) * 255.0) as u32;
                    (255 << 24) | (c << 16)
                }
                View::Charge => charge_to_color(self.charge[i], self.rho[i]),
                View::Matter | View::Potential => theta_rho_to_color(
                    self.v[i].atan2(self.u[i]),
                    self.rho[i],
                    matches!(view, View::Potential),
                ),
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f32, b: f32) {
        assert!((a - b).abs() < 1e-5, "{a} != {b}");
    }

    #[test]
    fn periodic_indices_wrap_in_both_directions() {
        let sim = Universe::new(7, 5);
        assert_eq!(sim.get_index(-1, -1), 34);
        assert_eq!(sim.get_index(7, 5), 0);
        assert_eq!(sim.get_index(-15, 11), sim.get_index(6, 1));
        for y in -10..10 {
            for x in -14..14 {
                assert_eq!(sim.get_index(x, y), sim.get_index(x + 21, y - 15));
            }
        }
    }

    #[test]
    fn uniform_field_stays_uniform_and_density_decays_by_existing_rule() {
        let mut sim = Universe::new(9, 7);
        sim.step();
        // Independent substitution: no gradients/stress; target_r=0, zero momentum.
        let expected_rho = 0.82 * (0.01 * 0.94) + 0.18 * 0.01;
        for i in 0..sim.u.len() {
            close(sim.u[i], 1.0);
            close(sim.v[i], 0.0);
            close(sim.rho[i], expected_rho);
            close(sim.charge[i], 0.0);
            close(sim.px[i], 0.0);
            close(sim.py[i], 0.0);
        }
        for _ in 0..100 {
            sim.step();
        }
        assert!(sim.is_finite());
        assert!(sim.u.iter().all(|&x| x == 1.0));
        assert!(sim.v.iter().all(|&x| x == 0.0));
        assert!(sim.rho.iter().all(|&x| x >= 0.001));
    }

    #[test]
    fn charge_measures_signed_winding_and_refreshes_after_step_and_spawn() {
        let mut sim = Universe::new(4, 4);
        let indices = [
            sim.get_index(1, 2),
            sim.get_index(2, 2),
            sim.get_index(2, 1),
            sim.get_index(1, 1),
        ];
        for (index, phase) in indices.into_iter().zip([0.0, PI / 2.0, PI, -PI / 2.0]) {
            sim.u[index] = phase.cos();
            sim.v[index] = phase.sin();
        }
        let charge = sim.recompute_charge_from_uv(&sim.u, &sim.v);
        close(charge[sim.get_index(1, 2)], 1.0);
        for v in &mut sim.v {
            *v = -*v;
        }
        let charge = sim.recompute_charge_from_uv(&sim.u, &sim.v);
        close(charge[sim.get_index(1, 2)], -1.0);
        let mut spawned = Universe::new(80, 80);
        spawned.spawn_vortex(40, 40, 1.0);
        assert_eq!(
            spawned.charge,
            spawned.recompute_charge_from_uv(&spawned.u, &spawned.v)
        );
        spawned.step();
        assert_eq!(
            spawned.charge,
            spawned.recompute_charge_from_uv(&spawned.u, &spawned.v)
        );
    }

    #[test]
    fn reset_restores_all_fields_and_history() {
        let mut sim = Universe::new(80, 64);
        sim.spawn_vortex(40, 32, -1.0);
        for _ in 0..10 {
            sim.step();
        }
        sim.push_nearby(30, 32, 0.8, 0.0);
        sim.reset();
        assert_eq!(sim, Universe::new(80, 64));
    }

    #[test]
    fn opposite_spawns_reverse_winding_in_current_and_previous_fields() {
        let mut positive = Universe::new(96, 80);
        let mut negative = positive.clone();
        positive.spawn_vortex(48, 40, 1.0);
        negative.spawn_vortex(48, 40, -1.0);
        // Both retain the existing PI/2 phase offset: u reverses, v remains.
        for i in 0..positive.u.len() {
            close(positive.u[i], -negative.u[i]);
            close(positive.v[i], negative.v[i]);
            close(positive.u_prev[i], -negative.u_prev[i]);
            close(positive.v_prev[i], negative.v_prev[i]);
        }
        // Core centers have zero amplitude, so inspect nonzero plaquettes adjacent to them.
        let positive_winding: f32 = (25..40)
            .flat_map(|x| (33..47).map(move |y| (x, y)))
            .map(|(x, y)| positive.charge[positive.get_index(x, y)])
            .sum();
        let negative_winding: f32 = (25..40)
            .flat_map(|x| (33..47).map(move |y| (x, y)))
            .map(|(x, y)| negative.charge[negative.get_index(x, y)])
            .sum();
        assert!(positive_winding.abs() > 0.9);
        close(positive_winding, -negative_winding);
    }

    #[test]
    fn deterministic_multistep_run_and_all_diagnostic_frames_remain_finite() {
        let mut sim = Universe::new(96, 80);
        sim.spawn_vortex(48, 40, 1.0);
        let mut repeat = sim.clone();
        for _ in 0..200 {
            sim.step();
            repeat.step();
            assert!(sim.is_finite());
        }
        assert_eq!(sim, repeat);
        for view in [
            View::Matter,
            View::Potential,
            View::Charge,
            View::Core,
            View::Momentum,
        ] {
            let pixels = sim.frame(view);
            assert_eq!(pixels.len(), 96 * 80);
            assert!(pixels.iter().all(|pixel| pixel >> 24 == 255));
        }
    }

    #[test]
    #[should_panic(expected = "grid dimensions must be positive")]
    fn empty_grid_is_rejected() {
        Universe::new(0, 5);
    }
}
