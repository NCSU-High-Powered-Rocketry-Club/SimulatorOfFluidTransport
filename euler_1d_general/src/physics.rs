use crate::GAMMA;
// ============================================================================
// 2. PHYSICS HELPER FUNCTIONS (Equation of State)
// ============================================================================

/// Calculates Pressure from Conserved Variables: p = (gamma - 1) * (E - 0.5 * rho * u^2)
/// State vector `u` is [rho, rho*u, E]
pub(crate) fn get_pressure(q: &[f64]) -> f64 {
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
pub(crate) fn get_sound_speed(q: &[f64]) -> f64 {
    let p = get_pressure(q);
    let rho = q[0];
    (GAMMA * p / rho).sqrt()
}



