//! # gooding-lambert
//!
//! Gooding's method for solving Lambert's orbital boundary-value problem.
//!
//! Given two position vectors and a time of flight, computes the velocity
//! vectors at each endpoint of the connecting Keplerian arc. Handles all
//! conic sections (elliptic, parabolic, hyperbolic), all transfer angles,
//! and multi-revolution solutions.
//!
//! ## Quick Start
//!
//! ```
//! use gooding_lambert::{lambert, Direction};
//!
//! let mu = 398600.4418_f64; // Earth GM (km³/s²)
//! let r1 = [6678.0, 0.0, 0.0]; // km
//! let r2 = [0.0, 42164.0, 0.0]; // km (GEO)
//! let tof = 5.0 * 3600.0; // seconds
//!
//! let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde).unwrap();
//! let speed = (sol.v1[0].powi(2) + sol.v1[1].powi(2) + sol.v1[2].powi(2)).sqrt();
//! assert!(speed > 5.0 && speed < 15.0); // km/s at LEO departure
//! ```
//!
//! ## References
//!
//! Gooding, R. H. (1990). "A procedure for the solution of Lambert's orbital
//! boundary-value problem." *Celestial Mechanics and Dynamical Astronomy*, 48(2), 145–165.

mod gooding;

pub use gooding::lambert;

/// Transfer direction for Lambert's problem.
///
/// Prograde means the transfer orbit has positive angular momentum aligned
/// with `r1 × r2`. Retrograde uses the supplementary transfer angle
/// (> 180°), with angular momentum aligned with `r2 × r1`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    /// Transfer angle ∈ (0°, 180°); angular momentum parallel to r1 × r2.
    Prograde,
    /// Transfer angle ∈ (180°, 360°); angular momentum parallel to r2 × r1.
    Retrograde,
}

/// Velocity solution returned by [`lambert`].
#[derive(Debug, Clone, PartialEq)]
pub struct LambertSolution {
    /// Velocity vector at the departure point (same units as position / time).
    pub v1: [f64; 3],
    /// Velocity vector at the arrival point (same units as position / time).
    pub v2: [f64; 3],
}

/// Errors returned by [`lambert`].
#[derive(Debug, Clone, PartialEq)]
pub enum LambertError {
    /// Transfer angle is exactly 180° — transfer plane undefined.
    SingularTransfer,
    /// No solution exists for the given revolution count and time of flight
    /// (TOF is below the minimum time for the requested revolution count).
    NoSolution,
    /// Householder iteration failed to converge within the iteration limit.
    ConvergenceFailed,
    /// One or more inputs are invalid (zero radius, non-positive TOF, etc.).
    InvalidInput(&'static str),
}
