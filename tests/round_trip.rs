//! Round-trip integration tests.
//!
//! Solve Lambert → propagate (r1, v1) forward by TOF → verify arrival at r2.

use gooding_lambert::{lambert, Direction};
use std::f64::consts::PI;

// ── Kepler propagator (universal variable + Stumpff) ───────────────────────

fn vec_mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn stumpff(psi: f64) -> (f64, f64) {
    if psi > 1e-6 {
        let sq = psi.sqrt();
        ((1.0 - sq.cos()) / psi, (sq - sq.sin()) / (psi * sq))
    } else if psi < -1e-6 {
        let sq = (-psi).sqrt();
        ((1.0 - sq.cosh()) / psi, (sq.sinh() - sq) / ((-psi) * sq))
    } else {
        (0.5 - psi / 24.0 + psi * psi / 720.0, 1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0)
    }
}

fn kepler_propagate(r0: [f64; 3], v0: [f64; 3], dt: f64, mu: f64) -> ([f64; 3], [f64; 3]) {
    let r0m = vec_mag(r0);
    let v0sq = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
    let rdotv = r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2];
    let alpha = 2.0 / r0m - v0sq / mu; // = 1/a (positive for elliptic)

    // Initial guess for universal variable chi
    let mut chi = if alpha > 1e-12 {
        mu.sqrt() * dt * alpha
    } else if alpha < -1e-12 {
        let a = 1.0 / alpha;
        let sign = if dt >= 0.0 { 1.0 } else { -1.0 };
        sign * (-a).sqrt()
            * ((-2.0 * mu * alpha * dt * dt)
                / (rdotv + sign * (-mu * a).sqrt() * (1.0 - r0m * alpha)))
                .ln()
    } else {
        mu.sqrt() * dt / r0m
    };

    // Newton–Raphson on Kepler's equation in universal variables
    let tol = 1e-14 * dt.abs().max(1.0);
    for _ in 0..50 {
        let psi = alpha * chi * chi;
        let (c2, c3) = stumpff(psi);
        let r = chi * chi * c2
            + rdotv / mu.sqrt() * chi * (1.0 - psi * c3)
            + r0m * (1.0 - psi * c2);
        let f_val = r0m * chi * (1.0 - psi * c3)
            + rdotv / mu.sqrt() * chi * chi * c2
            + chi * chi * chi * c3
            - mu.sqrt() * dt;
        let d = f_val / r;
        chi -= d;
        if d.abs() < tol {
            break;
        }
    }

    let psi = alpha * chi * chi;
    let (c2, c3) = stumpff(psi);
    let r_mag = chi * chi * c2
        + rdotv / mu.sqrt() * chi * (1.0 - psi * c3)
        + r0m * (1.0 - psi * c2);

    let f = 1.0 - chi * chi / r0m * c2;
    let g = dt - chi * chi * chi / mu.sqrt() * c3;
    let g_dot = 1.0 - chi * chi / r_mag * c2;
    let f_dot = mu.sqrt() / (r_mag * r0m) * chi * (psi * c3 - 1.0);

    let r = [
        f * r0[0] + g * v0[0],
        f * r0[1] + g * v0[1],
        f * r0[2] + g * v0[2],
    ];
    let v = [
        f_dot * r0[0] + g_dot * v0[0],
        f_dot * r0[1] + g_dot * v0[1],
        f_dot * r0[2] + g_dot * v0[2],
    ];
    (r, v)
}

// ── Test harness ────────────────────────────────────────────────────────────

fn round_trip(r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction, tol: f64) {
    let sol = lambert(mu, r1, r2, tof, 0, dir)
        .unwrap_or_else(|e| panic!("Lambert failed: {e:?}"));

    let (r2_prop, v2_prop) = kepler_propagate(r1, sol.v1, tof, mu);

    let pos_err = {
        let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
        vec_mag(d)
    };
    let r2_mag = vec_mag(r2);
    assert!(
        pos_err < tol * r2_mag,
        "position error {pos_err:.2e} (relative {:.2e}), tol={tol:.0e}",
        pos_err / r2_mag
    );

    // v2 agreement (propagated vs. Lambert)
    let vel_err = {
        let d = [
            sol.v2[0] - v2_prop[0],
            sol.v2[1] - v2_prop[1],
            sol.v2[2] - v2_prop[2],
        ];
        vec_mag(d)
    };
    let v2_mag = vec_mag(sol.v2);
    assert!(
        vel_err < tol * v2_mag,
        "velocity error {vel_err:.2e} (relative {:.2e}), tol={tol:.0e}",
        vel_err / v2_mag
    );
}

// ── Tests ───────────────────────────────────────────────────────────────────

#[test]
fn circular_90_deg() {
    round_trip([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn circular_45_deg() {
    let a = PI / 4.0;
    round_trip([1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn circular_150_deg() {
    let a = 150.0_f64.to_radians();
    round_trip([1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn elliptic_different_radii() {
    round_trip([1.0, 0.0, 0.0], [0.0, 1.5, 0.0], 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn hyperbolic() {
    round_trip([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn retrograde() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        3.0 * PI / 2.0,
        1.0,
        Direction::Retrograde,
        1e-8,
    );
}

#[test]
fn three_dimensional() {
    round_trip([1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0, Direction::Prograde, 1e-8);
}

#[test]
fn physical_units_leo_to_geo() {
    let mu = 398600.4418;
    round_trip(
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        5.0 * 3600.0,
        mu,
        Direction::Prograde,
        1e-8,
    );
}
