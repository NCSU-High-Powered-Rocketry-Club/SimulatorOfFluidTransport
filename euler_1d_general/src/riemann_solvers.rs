use ndarray::{Array1,array};

use crate::constants::*;
use crate::physics::{get_sound_speed,get_pressure};

pub(crate) fn ldfss_flux(uvec_l: &[f64], uvec_r: &[f64]) -> Array1::<f64> {
    // 1. Primitive Variables
    let rho_l = uvec_l[0];
    let u_l   = uvec_l[1] / rho_l;
    let p_l   = get_pressure(uvec_l);
    let h_l   = (uvec_l[2] + p_l) / rho_l; // Total Enthalpy
    let a_l   = get_sound_speed(uvec_l); // Sound speed

    let rho_r = uvec_r[0];
    let u_r   = uvec_r[1] / rho_r;
    let p_r   = get_pressure(uvec_r);
    let h_r   = (uvec_r[2] + p_r) / rho_r;
    let a_r   = get_sound_speed(uvec_r); // Sound speed

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
