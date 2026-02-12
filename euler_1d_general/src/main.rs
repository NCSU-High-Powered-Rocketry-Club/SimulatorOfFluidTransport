use ndarray::{Array1,Array2,array};
use std::fs::File;
use std::io::Write;

// ============================================================================
// 1. PHYSICAL CONSTANTS & CONFIGURATION
// ============================================================================
const GAMMA: f64 = 1.4;          // Ratio of specific heats for ideal gas
const N_CELLS: usize = 400;      // Number of spatial grid points
const L_DOMAIN: f64 = 1.0;       // Length of the domain (meters)
const T_FINAL: f64 = 0.2;        // Final simulation time (seconds)
const CFL: f64 = 0.5;            // Courant-Friedrichs-Lewy stability number

// ============================================================================
// 2. PHYSICS HELPER FUNCTIONS (Equation of State)
// ============================================================================

/// Calculates Pressure from Conserved Variables: p = (gamma - 1) * (E - 0.5 * rho * u^2)
/// State vector `u` is [rho, rho*u, E]
fn get_pressure(q: &[f64]) -> f64 {
    let rho = q[0];
    let mom = q[1];
    let energy = q[2];
    
    // Safety check for vacuum or invalid state
    if rho <= 1.0e-10 { return 1.0e-10; }
    
    let u_velocity = mom / rho;
    let kinetic_energy = 0.5 * rho * u_velocity * u_velocity;
    let internal_energy = energy - kinetic_energy;
    
    (GAMMA - 1.0) * internal_energy
}

/// Calculates Sound Speed: c = sqrt(gamma * p / rho)
fn get_sound_speed(q: &[f64]) -> f64 {
    let p = get_pressure(q);
    let rho = q[0];
    (GAMMA * p / rho).sqrt()
}

// ============================================================================
// 3. NUMERICAL FLUX: EDWARDS' LDFSS (Low-Diffusion Flux-Splitting Scheme)
// ============================================================================
// Note: This implements the core AUSM+ logic typical of Edwards' LDFSS.
// It splits fluxes into Convective (Mass/Velocity) and Pressure parts separately.

fn ldfss_flux(q_l: &[f64], q_r: &[f64]) -> Array1::<f64> {
    // 1. Primitive Variables
    let rho_l = q_l[0];
    let u_l   = q_l[1] / rho_l;
    let p_l   = get_pressure(q_l);
    let h_l   = (q_l[2] + p_l) / rho_l; // Total Enthalpy
    let a_l   = (GAMMA * p_l / rho_l).sqrt(); // Sound speed

    let rho_r = q_r[0];
    let u_r   = q_r[1] / rho_r;
    let p_r   = get_pressure(q_r);
    let h_r   = (q_r[2] + p_r) / rho_r;
    let a_r   = (GAMMA * p_r / rho_r).sqrt();

    // 2. Interface Sound Speed (Simple arithmetic mean often used in LDFSS)
    let a_half = 0.5 * (a_l + a_r);

    // 3. Mach Numbers at Interface
    let m_l = u_l / a_half;
    let m_r = u_r / a_half;

    // 4. Split Mach Polynomials M+ (Left) and M- (Right)
    // Using standard Van Leer / AUSM polynomials
    let m_plus = if m_l.abs() >= 1.0 {
        0.5 * (m_l + m_l.abs())
    } else {
        0.25 * (m_l + 1.0).powi(2)
    };

    let m_minus = if m_r.abs() >= 1.0 {
        0.5 * (m_r - m_r.abs())
    } else {
        -0.25 * (m_r - 1.0).powi(2)
    };

    // 5. Interface Mach Number
    let m_half = m_plus + m_minus;

    // 6. Split Pressure Polynomials P+ and P-
    let p_plus = if m_l.abs() >= 1.0 {
        if m_l > 0.0 { 1.0 } else { 0.0 }
    } else {
        0.25 * (m_l + 1.0).powi(2) * (2.0 - m_l)
    };

    let p_minus = if m_r.abs() >= 1.0 {
        if m_r > 0.0 { 0.0 } else { 1.0 }
    } else {
        0.25 * (m_r - 1.0).powi(2) * (2.0 + m_r)
    };

    // 7. Assemble Fluxes
    // F = F_convective + F_pressure
    // The convective flux is determined by the sign of the Interface Mach Number (Upwinding)
    
    let mass_flux = a_half * m_half; 

    // Compute convective vector based on upwinding direction of mass_flux
    let (rho_up, u_up, h_up) = if mass_flux >= 0.0 {
        (rho_l, u_l, h_l)
    } else {
        (rho_r, u_r, h_r)
    };

    // Convective Flux Terms: [rho*u, rho*u^2, rho*u*H]
    let f_conv_0 = mass_flux * rho_up;
    let f_conv_1 = mass_flux * rho_up * u_up;
    let f_conv_2 = mass_flux * rho_up * h_up;

    // Pressure Flux Terms: [0, p, 0]
    // The pressure at interface is weighted blend of P_L and P_R
    let p_half = p_plus * p_l + p_minus * p_r;

    // Total Flux
    array![
        f_conv_0,
        f_conv_1 + p_half,
        f_conv_2 // Enthalpy flux contains the work term
    ]
}

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
