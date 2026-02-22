use ndarray::{Array1,Array2,ArrayViewMut1};
use std::fs::File;
use std::io::Write;


mod constants;
mod riemann_solvers;
mod physics;
mod type_solution;

use crate::constants::*;
use crate::riemann_solvers::riemann_solver;
use crate::physics::{get_sound_speed,get_pressure};
use crate::type_solution::SolutionArray;

// ============================================================================
// MAIN SOLVER
// ============================================================================

fn main() -> std::io::Result<()> {
    // --- A. Grid Initialization ---
    let dx = L_DOMAIN / (N_CELLS as f64);
    // 
    // We use Array2 for the state vector U. 
    // Shape: (N_CELLS, 3) -> Rows are cells, Columns are [rho, mom, E]
    let mut u      = SolutionArray::new(N_CELLS,   NVARS);
    let mut u_new  = SolutionArray::new(N_CELLS,   NVARS);
    let mut fluxes = SolutionArray::new(N_CELLS+1, NVARS);
    
    // --- B. Initial Conditions (Sod Shock Tube) ---
    for i in 0..N_CELLS {
        let x = (i as f64 + 0.5) * dx;
        initialize_sod_shock(x, u.row_mut(i))
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
        for i in 1..N_CELLS {
            // Get Left and Right states for the interface between i and i+1
            let u_l = u.row(i-1);
            let u_r = u.row(i);
            
            // Compute the common interface flux values
            fluxes.row_mut(i).assign(
                &riemann_solver(
                    u_l.as_slice().unwrap(), 
                    u_r.as_slice().unwrap(),
                    RIEMANN_OPT)
            );
        }
        // Boundary Conditions (Transmissive / Zero Gradient)
        // Flux at 0 (Left Boundary) computed using cell 0 as both L and R (approximated)
        // or just copying flux from neighbor. 
        fluxes.row_mut(0).assign( 
            &riemann_solver(
                u.row(0).as_slice().unwrap(), 
                u.row(0).as_slice().unwrap(),
                RIEMANN_OPT)
        );
        fluxes.row_mut(N_CELLS).assign(
            &riemann_solver(
                u.row(N_CELLS-1).as_slice().unwrap(),
                u.row(N_CELLS-1).as_slice().unwrap(),
                RIEMANN_OPT)
        );

        // 3. Update Solution (Finite Volume Formulation)
        // U_new = U - dt/dx * (F_right - F_left)
        for i in 2..N_CELLS-2 {
            let f_right = fluxes.row(i+1);
            let f_left  = fluxes.row(i);
            u_new.row_mut(i).assign(  &(
                &u.row(i) - (dt/dx)*(&f_right - &f_left)
            ));
        }
        
        // hardcode pinned BC states (Zero Gradient)
        let mut u_ghost = u_new.row(1).to_owned();
        u_new.row_mut(0).assign(&u_ghost);
        
        u_ghost = u_new.row(N_CELLS-2).to_owned();
        u_new.row_mut(N_CELLS-1).assign(&u_ghost);

        // Swap buffers
        u.data.assign(&u_new.data);
        
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
        let rho    = u.data[[i, NDENS]];
        let mom    = u.data[[i, NMOMX]];
        let energy = u.data[[i, NERGY]];
        
        let velocity = mom / rho;
        let p = get_pressure(u.row(i).as_slice().unwrap());

        writeln!(file, "{:.6},{:.6},{:.6},{:.6},{:.6}", x, rho, velocity, p, energy)?;
    }

    Ok(())
}
//
pub (crate) fn initialize_sod_shock(x: f64, mut uout: ArrayViewMut1<f64>) {
    // Left:  Rho=1.0,   P=1.0, u=0.0
    // Right: Rho=0.125, P=0.1, u=0.0
    // Diaphragm at x = 0.5
    
    let (rho, p, velocity) = if x < 0.5 * L_DOMAIN {
        (1.0, 1.0, 0.0)
    } else {
        (0.125, 0.1, 0.0)
    };

    let mom = rho * velocity;
    let energy = p / (GAMMA - 1.0) + 0.5 * rho * velocity * velocity;

    uout[NDENS] = rho;
    uout[NMOMX] = mom;
    uout[NERGY] = energy;
}
