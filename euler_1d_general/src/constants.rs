// ============================================================================
// 1. PHYSICAL CONSTANTS & CONFIGURATION
// ============================================================================
pub const GAMMA: f64 = 1.4;          // Ratio of specific heats for ideal gas
pub const N_VARIABLES_EULER_1D: usize = 3;
//
// Simulation Controls (should move to inputs eventually)
pub const N_CELLS: usize = 400;      // Number of spatial grid points
pub const L_DOMAIN: f64 = 1.0;       // Length of the domain (meters)
pub const T_FINAL: f64 = 0.2;        // Final simulation time (seconds)
pub const CFL: f64 = 0.5;            // Courant-Friedrichs-Lewy stability number
