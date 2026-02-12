// ============================================================================
// 1. PHYSICAL CONSTANTS & CONFIGURATION
// ============================================================================
pub const GAMMA: f64 = 1.4;          // Ratio of specific heats for ideal gas
pub const N_CELLS: usize = 400;      // Number of spatial grid points
pub const L_DOMAIN: f64 = 1.0;       // Length of the domain (meters)
pub const T_FINAL: f64 = 0.2;        // Final simulation time (seconds)
pub const CFL: f64 = 0.5;            // Courant-Friedrichs-Lewy stability number
//
pub const AS_IMM: i32 = 0;
pub const AS_MUT: i32 = 1;
pub const AS_OWN: i32 = 2;
