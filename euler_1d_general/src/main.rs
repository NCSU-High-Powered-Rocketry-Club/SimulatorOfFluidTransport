use ndarray::{Array1,Array2};
use std::fs::File;
use std::io::Write;


mod constants;
mod riemann_solvers;
mod physics;

use crate::constants::*;
use crate::riemann_solvers::ldfss_flux;
use crate::physics::{get_sound_speed,get_pressure};




// ============================================================================
// 4. MAIN SOLVER
// ============================================================================

fn main() -> std::io::Result<()> {
    // --- A. Grid Initialization ---
    let dx = L_DOMAIN / (N_CELLS as f64);
    
    // We use Array2 for the state vector U. 
    // Shape: (N_CELLS, 3) -> Rows are cells, Columns are [rho, mom, E]
    let mut u      = Array2::<f64>::zeros((N_CELLS, 3));
    let mut u_new  = Array2::<f64>::zeros((N_CELLS, 3));
    let mut fluxes = Array2::<f64>::zeros((N_CELLS+1, 3));
    
    // --- B. Initial Conditions (Sod Shock Tube) ---
    // Left:  Rho=1.0,   P=1.0, u=0.0
    // Right: Rho=0.125, P=0.1, u=0.0
    // Diaphragm at x = 0.5
    for i in 0..N_CELLS {
        let x = (i as f64 + 0.5) * dx;
        
        let (rho, p, velocity) = if x < 0.5 * L_DOMAIN {
            (1.0, 1.0, 0.0)
        } else {
            (0.125, 0.1, 0.0)
        };

        let mom = rho * velocity;
        let energy = p / (GAMMA - 1.0) + 0.5 * rho * velocity * velocity;

        u[[i, 0]] = rho;
        u[[i, 1]] = mom;
        u[[i, 2]] = energy;
    }

    // --- C. Time Loop ---
    let mut t = 0.0;
    let mut step = 0;

    println!("Starting Simulation...");
    println!("Grid: {}, dx: {:.6}, CFL: {}", N_CELLS, dx, CFL);

    while t < T_FINAL {
        // 1. Calculate Time Step (dt) based on max wave speed
        let mut max_speed = 0.0;
        for i in 0..N_CELLS {
            let q = u.row(i);
            let rho = q[0];
            let u_vel = (q[1] / rho).abs();
            let c = get_sound_speed(q.as_slice().unwrap());
            let wave_speed = u_vel + c;
            if wave_speed > max_speed {
                max_speed = wave_speed;
            }
        }
        
        let mut dt = CFL * dx / max_speed;
        if t + dt > T_FINAL {
            dt = T_FINAL - t;
        }

        // 2. Compute Fluxes at Interfaces
        // We have N cells, so we have N+1 interfaces.
        // Interface i is between Cell i-1 and Cell i.
        // We store fluxes in a generic Vec for simplicity or an Array2.
        // Let's use a Vec of [f64; 3] to hold fluxes at i+1/2.
        // Indices: 0 to N. Index i corresponds to interface i-1/2 ?? 
        // Let's stick to: flux[i] is flux at right face of cell i.
        

        for i in 0..N_CELLS-1 {
            // Get Left and Right states for the interface between i and i+1
            let u_l = u.row(i);
            let u_r = u.row(i + 1);
            
            fluxes.row_mut(i+1)  .assign( &ldfss_flux(u_l.as_slice().unwrap(), u_r.as_slice().unwrap())  );
        }

        // Boundary Conditions (Transmissive / Zero Gradient)
        // Flux at 0 (Left Boundary) computed using cell 0 as both L and R (approximated)
        // or just copying flux from neighbor. 
        // Standard Sod: Left/Right boundaries don't interact within T=0.2, 
        // so we can just impose static fluxes or zero gradient.
        fluxes.row_mut(0)       .assign( &ldfss_flux(u.row(0).as_slice().unwrap(), u.row(0).as_slice().unwrap()) );
        fluxes.row_mut(N_CELLS) .assign( &ldfss_flux(u.row(N_CELLS-1).as_slice().unwrap(), u.row(N_CELLS-1).as_slice().unwrap()) );

        // 3. Update Solution (Finite Volume Formulation)
        // U_new = U - dt/dx * (F_right - F_left)
        for i in 1..N_CELLS-1 {
            let f_right = fluxes.row(i+1);
            let f_left  = fluxes.row(i);
            
            for var in 0..3 {
                u_new[[i, var]] = u[[i, var]] - (dt / dx) * (f_right[var] - f_left[var]);
            }
        }

        // Apply BCs to states (Zero Gradient)
        for var in 0..3 {
            u_new[[0, var]] = u_new[[1, var]];
            u_new[[N_CELLS-1, var]] = u_new[[N_CELLS-2, var]];
        }

        // Swap buffers
        u.assign(&u_new);
        
        t += dt;
        step += 1;
        
        if step % 50 == 0 {
            println!("Step: {}, Time: {:.4}, dt: {:.6}", step, t, dt);
        }
    }

    // --- D. Output Results ---
    println!("Simulation Complete. Writing results to 'sod_results.csv'");
    let mut file = File::create("sod_results.csv")?;
    writeln!(file, "x,density,velocity,pressure,energy")?;

    for i in 0..N_CELLS {
        let x = (i as f64 + 0.5) * dx;
        let rho = u[[i, 0]];
        let mom = u[[i, 1]];
        let energy = u[[i, 2]];
        
        let velocity = mom / rho;
        let p = get_pressure(u.row(i).as_slice().unwrap());

        writeln!(file, "{:.6},{:.6},{:.6},{:.6},{:.6}", x, rho, velocity, p, energy)?;
    }

    Ok(())
}
