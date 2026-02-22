// ============================================================================
// 1. PHYSICAL CONSTANTS & CONFIGURATION
// ============================================================================
pub const GAMMA: f64 = 1.4;          // Ratio of specific heats for ideal gas
pub const N_VARIABLES_EULER_1D: usize = 3;
// Riemann / Interface flux options
pub const LDFSS: usize = 0;
pub const AUSM: usize = 1;
//
// Simulation Controls (should move to inputs eventually)
pub const N_CELLS: usize = 400;      // Number of spatial grid points
pub const L_DOMAIN: f64 = 1.0;       // Length of the domain (meters)
pub const T_FINAL: f64 = 0.2;        // Final simulation time (seconds)
pub const CFL: f64 = 0.5;            // Courant-Friedrichs-Lewy stability number
// --- Case specific numbers setup ---
pub const NVARS: usize = N_VARIABLES_EULER_1D;
pub const NDENS: usize = 0;
pub const NMOMX: usize = 1;
pub const NERGY: usize = 2;
pub const RIEMANN_OPT: usize = LDFSS;
