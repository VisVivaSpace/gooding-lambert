//! Round-trip integration tests.
//!
//! Solve Lambert -> propagate (r1, v1) forward by TOF -> verify arrival at r2.

mod common;

use gooding_lambert::{Direction, MultiRevPeriod, lambert};
use std::f64::consts::PI;

// -- Kepler propagator (universal variable + Stumpff) -------------------------

fn stumpff(psi: f64) -> (f64, f64) {
    if psi > 1e-6 {
        let sq = psi.sqrt();
        ((1.0 - sq.cos()) / psi, (sq - sq.sin()) / (psi * sq))
    } else if psi < -1e-6 {
        let sq = (-psi).sqrt();
        ((1.0 - sq.cosh()) / psi, (sq.sinh() - sq) / ((-psi) * sq))
    } else {
        (
            0.5 - psi / 24.0 + psi * psi / 720.0,
            1.0 / 6.0 - psi / 120.0 + psi * psi / 5040.0,
        )
    }
}

/// Kepler propagation via universal variables.
///
/// Returns `Some((r, v))` on convergence, `None` if Newton-Raphson diverges.
fn kepler_propagate(r0: [f64; 3], v0: [f64; 3], dt: f64, mu: f64) -> Option<([f64; 3], [f64; 3])> {
    let r0m = common::vec_mag(r0);
    let v0sq = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
    let rdotv = r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2];
    let alpha = 2.0 / r0m - v0sq / mu; // = 1/a (positive for elliptic)

    // Initial guess for universal variable chi
    let mut chi = if alpha > 1e-12 {
        mu.sqrt() * dt * alpha
    } else if alpha < -1e-12 {
        let a = 1.0 / alpha;
        let sign = if dt >= 0.0 { 1.0 } else { -1.0 };
        let arg = (-2.0 * mu * alpha * dt * dt)
            / (rdotv + sign * (-mu * a).sqrt() * (1.0 - r0m * alpha));
        if !arg.is_finite() || arg <= 0.0 {
            return None;
        }
        sign * (-a).sqrt() * arg.ln()
    } else {
        mu.sqrt() * dt / r0m
    };

    if !chi.is_finite() {
        return None;
    }

    // Newton-Raphson on Kepler's equation in universal variables
    let tol = 1e-14 * dt.abs().max(1.0);
    let mut converged = false;
    for _ in 0..50 {
        let psi = alpha * chi * chi;
        let (c2, c3) = stumpff(psi);
        let r =
            chi * chi * c2 + rdotv / mu.sqrt() * chi * (1.0 - psi * c3) + r0m * (1.0 - psi * c2);
        if !r.is_finite() || r.abs() < 1e-30 {
            return None;
        }
        let f_val = r0m * chi * (1.0 - psi * c3)
            + rdotv / mu.sqrt() * chi * chi * c2
            + chi * chi * chi * c3
            - mu.sqrt() * dt;
        let d = f_val / r;
        chi -= d;
        if !chi.is_finite() {
            return None;
        }
        if d.abs() < tol {
            converged = true;
            break;
        }
    }

    if !converged {
        return None;
    }

    let psi = alpha * chi * chi;
    let (c2, c3) = stumpff(psi);
    let r_mag =
        chi * chi * c2 + rdotv / mu.sqrt() * chi * (1.0 - psi * c3) + r0m * (1.0 - psi * c2);

    if !r_mag.is_finite() || r_mag.abs() < 1e-30 {
        return None;
    }

    let f = 1.0 - chi * chi / r0m * c2;
    let g = dt - chi * chi * chi / mu.sqrt() * c3;
    let g_dot = 1.0 - chi * chi / r_mag * c2;
    let f_dot = mu.sqrt() / (r_mag * r0m) * chi * (psi * c3 - 1.0);

    // Verify Lagrange identity: f*g_dot - f_dot*g = 1
    let identity = f * g_dot - f_dot * g;
    if (identity - 1.0).abs() > 1e-8 {
        return None;
    }

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
    Some((r, v))
}

// -- Test harness -------------------------------------------------------------

fn round_trip(r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction, tol: f64) {
    let sol = lambert(mu, r1, r2, tof, 0, dir, MultiRevPeriod::LongPeriod).unwrap_or_else(|e| panic!("Lambert failed: {e:?}"));

    let (r2_prop, v2_prop) = kepler_propagate(r1, sol.v1, tof, mu)
        .expect("Kepler propagator failed to converge");

    let pos_err = {
        let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
        common::vec_mag(d)
    };
    let r2_mag = common::vec_mag(r2);
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
        common::vec_mag(d)
    };
    let v2_mag = common::vec_mag(sol.v2);
    assert!(
        vel_err < tol * v2_mag,
        "velocity error {vel_err:.2e} (relative {:.2e}), tol={tol:.0e}",
        vel_err / v2_mag
    );
}

// -- Tests --------------------------------------------------------------------
//
// Tolerance derivation (Tier 3 -- iterative method):
//   xlamb converges to |dT| < 1e-12 (relative to T).
//   x -> velocity through ~15 arithmetic ops: rounding error ~15e ~= 3e-15.
//   Kepler propagator converges to |dchi| < 1e-14 * dt.max(1) (absolute).
//   For the shortest test (dt = 0.1), that is 1e-14; for LEO->GEO (dt = 18 000),
//   the chi error translates to position error ~3e-12 relative (see analysis in
//   docs/algorithms/gooding-lambert.md).
//   Net expected round-trip error: ~1e-12 relative.
//   Test tolerance: 1e-11 (10x safety margin over expected 1e-12).

#[test]
fn circular_90_deg() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        PI / 2.0,
        1.0,
        Direction::Prograde,
        1e-11,
    );
}

#[test]
fn circular_45_deg() {
    let a = PI / 4.0;
    round_trip(
        [1.0, 0.0, 0.0],
        [a.cos(), a.sin(), 0.0],
        a,
        1.0,
        Direction::Prograde,
        1e-11,
    );
}

#[test]
fn circular_150_deg() {
    let a = 150.0_f64.to_radians();
    round_trip(
        [1.0, 0.0, 0.0],
        [a.cos(), a.sin(), 0.0],
        a,
        1.0,
        Direction::Prograde,
        1e-11,
    );
}

#[test]
fn elliptic_different_radii() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        2.0,
        1.0,
        Direction::Prograde,
        1e-11,
    );
}

#[test]
fn hyperbolic() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        0.1,
        1.0,
        Direction::Prograde,
        1e-11,
    );
}

#[test]
fn retrograde() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        PI / 2.0,
        1.0,
        Direction::Retrograde,
        1e-11,
    );
}

#[test]
fn three_dimensional() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 0.8, 0.6],
        PI / 2.0,
        1.0,
        Direction::Prograde,
        1e-11,
    );
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
        1e-11,
    );
}

// -- Retrograde round-trip tests ----------------------------------------------

#[test]
fn retrograde_3d() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 0.8, 0.6],
        PI / 2.0,
        1.0,
        Direction::Retrograde,
        1e-11,
    );
}

#[test]
fn retrograde_different_radii() {
    round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.5, 0.0],
        2.0,
        1.0,
        Direction::Retrograde,
        1e-11,
    );
}

#[test]
fn retrograde_135_deg() {
    let a = 135.0_f64.to_radians();
    round_trip(
        [1.0, 0.0, 0.0],
        [a.cos(), a.sin(), 0.0],
        a,
        1.0,
        Direction::Retrograde,
        1e-11,
    );
}

#[test]
fn retrograde_leo_to_geo() {
    let mu = 398600.4418;
    round_trip(
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        5.0 * 3600.0,
        mu,
        Direction::Retrograde,
        1e-11,
    );
}

// -- Randomized round-trip tests ----------------------------------------------
//
// Generate random geometries, solve Lambert, propagate with Kepler, verify
// arrival at r2.  These tests check physical correctness -- a bug that
// produces a self-consistent wrong answer in both v1 and v2 would still fail
// here because the Kepler propagator is an independent implementation.
//
// Seeds are truly random (time-based).  On failure the seed is printed via
// eprintln! so the case can be reproduced with `cargo test -- --nocapture`.

/// Run `n` randomized round-trip cases for a given direction and mu.
///
/// Returns the number of cases actually tested (some are skipped when the
/// solver returns `NoSolution` for extreme geometries).
fn randomized_round_trips(
    seed: &mut u64,
    n: usize,
    mu: f64,
    r_lo: f64,
    r_hi: f64,
    dir: Direction,
    tol: f64,
    label: &str,
) -> usize {
    const NEAR_COLLINEAR: f64 = 0.9962; // cos(5 deg)

    let mut tested = 0;
    let mut attempts = 0;

    while tested < n {
        attempts += 1;
        assert!(
            attempts < n * 100,
            "{label}: gave up after {attempts} attempts (only {tested} tested); seed={seed}"
        );

        let r1 = common::rand_r(seed, r_lo, r_hi);
        let r2 = common::rand_r(seed, r_lo, r_hi);

        // Reject near-collinear geometries (transfer angle near 0 deg or 180 deg)
        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        // TOF: fraction of the *minimum-energy* orbital period.
        // Use the larger radius to set scale -- this avoids near-parabolic
        // solutions that would stress the Kepler propagator beyond its
        // convergence regime.
        let r1_mag = common::vec_mag(r1);
        let r2_mag = common::vec_mag(r2);
        let r_max = r1_mag.max(r2_mag);
        let t_circ = 2.0 * PI * (r_max * r_max * r_max / mu).sqrt();
        let tof = common::rand_f64(seed, 0.05 * t_circ, 0.8 * t_circ);

        // Solve Lambert -- skip cases where no solution exists
        let sol = match lambert(mu, r1, r2, tof, 0, dir, MultiRevPeriod::LongPeriod) {
            Ok(s) => s,
            Err(_) => continue,
        };

        // Propagate (r1, v1) forward by tof using Kepler -- skip if
        // the propagator does not converge (near-parabolic edge cases).
        let (r2_prop, v2_prop) = match kepler_propagate(r1, sol.v1, tof, mu) {
            Some(rv) => rv,
            None => continue,
        };

        // Check position agreement
        let pos_err = {
            let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
            common::vec_mag(d)
        };
        let r2_mag = common::vec_mag(r2);
        assert!(
            pos_err < tol * r2_mag,
            "{label}: position error {pos_err:.2e} (relative {:.2e}), tol={tol:.0e}\n  \
             seed={seed}, attempt={attempts}\n  \
             r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}\n  v1={:.6?}",
            pos_err / r2_mag,
            sol.v1
        );

        // Check velocity agreement
        let vel_err = {
            let d = [
                sol.v2[0] - v2_prop[0],
                sol.v2[1] - v2_prop[1],
                sol.v2[2] - v2_prop[2],
            ];
            common::vec_mag(d)
        };
        let v2_mag = common::vec_mag(sol.v2);
        assert!(
            vel_err < tol * v2_mag,
            "{label}: velocity error {vel_err:.2e} (relative {:.2e}), tol={tol:.0e}\n  \
             seed={seed}, attempt={attempts}\n  \
             r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}\n  v1={:.6?}",
            vel_err / v2_mag,
            sol.v1
        );

        tested += 1;
    }

    tested
}

#[test]
fn random_prograde_100() {
    let mut seed = common::make_seed();
    eprintln!("random_prograde_100: seed = {seed}");

    let tested = randomized_round_trips(
        &mut seed,
        100,
        1.0,   // canonical mu
        0.5,   // r_lo
        5.0,   // r_hi
        Direction::Prograde,
        1e-11,
        "random_prograde_100",
    );

    eprintln!("random_prograde_100: {tested} cases passed (seed={seed})");
}

#[test]
fn random_retrograde_100() {
    let mut seed = common::make_seed();
    eprintln!("random_retrograde_100: seed = {seed}");

    let tested = randomized_round_trips(
        &mut seed,
        100,
        1.0,   // canonical mu
        0.5,   // r_lo
        5.0,   // r_hi
        Direction::Retrograde,
        1e-11,
        "random_retrograde_100",
    );

    eprintln!("random_retrograde_100: {tested} cases passed (seed={seed})");
}

#[test]
fn random_physical_units_100() {
    let mut seed = common::make_seed();
    eprintln!("random_physical_units_100: seed = {seed}");

    let mu = 398600.4418; // Earth GM (km^3/s^2)

    // Mix prograde and retrograde in physical units (LEO to GEO range)
    let prograde = randomized_round_trips(
        &mut seed,
        50,
        mu,
        6400.0,   // r_lo: just above Earth radius (km)
        50000.0,  // r_hi: beyond GEO (km)
        Direction::Prograde,
        1e-11,
        "random_physical_units_prograde",
    );
    let retrograde = randomized_round_trips(
        &mut seed,
        50,
        mu,
        6400.0,
        50000.0,
        Direction::Retrograde,
        1e-11,
        "random_physical_units_retrograde",
    );

    eprintln!(
        "random_physical_units_100: {prograde} prograde + {retrograde} retrograde passed (seed={seed})"
    );
}

#[test]
fn random_high_eccentricity_50() {
    let mut seed = common::make_seed();
    eprintln!("random_high_eccentricity_50: seed = {seed}");

    // Wide radius range produces high-eccentricity transfers
    let tested = randomized_round_trips(
        &mut seed,
        50,
        1.0,
        0.2,    // r_lo: small periapsis
        10.0,   // r_hi: large apoapsis => eccentric transfers
        Direction::Prograde,
        1e-11,
        "random_high_eccentricity_50",
    );

    eprintln!("random_high_eccentricity_50: {tested} cases passed (seed={seed})");
}

#[test]
fn random_3d_out_of_plane_50() {
    let mut seed = common::make_seed();
    eprintln!("random_3d_out_of_plane_50: seed = {seed}");

    // rand_r already generates full 3D random directions -- this test
    // specifically verifies that out-of-plane geometries (large z-components)
    // work correctly.  We generate positions and filter for cases where both
    // vectors have significant z-components (|z| > 0.3 * |r|).
    const NEAR_COLLINEAR: f64 = 0.9962;
    let mut tested = 0;
    let mut attempts = 0;

    while tested < 50 {
        attempts += 1;
        assert!(
            attempts < 5000,
            "random_3d_out_of_plane_50: gave up after {attempts} attempts; seed={seed}"
        );

        let r1 = common::rand_r(&mut seed, 0.5, 5.0);
        let r2 = common::rand_r(&mut seed, 0.5, 5.0);

        // Require significant out-of-plane component
        let r1_mag = common::vec_mag(r1);
        let r2_mag = common::vec_mag(r2);
        if r1[2].abs() < 0.3 * r1_mag || r2[2].abs() < 0.3 * r2_mag {
            continue;
        }
        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        let r_mean = (r1_mag + r2_mag) / 2.0;
        let t_circ = 2.0 * PI * (r_mean * r_mean * r_mean / 1.0_f64).sqrt();
        let tof = common::rand_f64(&mut seed, 0.05 * t_circ, 1.5 * t_circ);

        let dir = if common::xorshift(&mut seed) % 2 == 0 {
            Direction::Prograde
        } else {
            Direction::Retrograde
        };

        let sol = match lambert(1.0, r1, r2, tof, 0, dir, MultiRevPeriod::LongPeriod) {
            Ok(s) => s,
            Err(_) => continue,
        };

        let (r2_prop, v2_prop) = match kepler_propagate(r1, sol.v1, tof, 1.0) {
            Some(rv) => rv,
            None => continue,
        };

        let pos_err = {
            let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
            common::vec_mag(d)
        };
        assert!(
            pos_err < 1e-11 * r2_mag,
            "random_3d_out_of_plane_50: position error {pos_err:.2e} (relative {:.2e})\n  \
             seed={seed}, attempt={attempts}\n  r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}",
            pos_err / r2_mag
        );

        let vel_err = {
            let d = [
                sol.v2[0] - v2_prop[0],
                sol.v2[1] - v2_prop[1],
                sol.v2[2] - v2_prop[2],
            ];
            common::vec_mag(d)
        };
        let v2_mag = common::vec_mag(sol.v2);
        assert!(
            vel_err < 1e-11 * v2_mag,
            "random_3d_out_of_plane_50: velocity error {vel_err:.2e} (relative {:.2e})\n  \
             seed={seed}, attempt={attempts}\n  r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}",
            vel_err / v2_mag
        );

        tested += 1;
    }

    eprintln!("random_3d_out_of_plane_50: {tested} cases passed in {attempts} attempts (seed={seed})");
}

// -- Multi-revolution round-trip tests ----------------------------------------
//
// Solve Lambert for nrev >= 1, propagate (r1, v1) forward by TOF with
// Kepler, verify arrival at r2.  Tests both LongPeriod and ShortPeriod
// solution families.
//
// For multi-rev solutions to exist, the TOF must exceed a minimum time
// that depends on the geometry and revolution count.  We use generous
// TOFs (several orbital periods) to stay well above the minimum.
//
// Tolerance: multi-rev solutions traverse more revolutions and accumulate
// more Kepler-propagation rounding, so we relax to 1e-9 (vs 1e-11 for
// zero-rev).

fn multi_rev_round_trip(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    nrev: u32,
    dir: Direction,
    period: MultiRevPeriod,
    tol: f64,
    label: &str,
) {
    let sol = lambert(mu, r1, r2, tof, nrev, dir, period)
        .unwrap_or_else(|e| panic!("{label}: Lambert failed: {e:?}"));

    let (r2_prop, v2_prop) = kepler_propagate(r1, sol.v1, tof, mu)
        .unwrap_or_else(|| panic!("{label}: Kepler propagator failed to converge"));

    let pos_err = {
        let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
        common::vec_mag(d)
    };
    let r2_mag = common::vec_mag(r2);
    assert!(
        pos_err < tol * r2_mag,
        "{label}: position error {pos_err:.2e} (relative {:.2e}), tol={tol:.0e}",
        pos_err / r2_mag
    );

    let vel_err = {
        let d = [
            sol.v2[0] - v2_prop[0],
            sol.v2[1] - v2_prop[1],
            sol.v2[2] - v2_prop[2],
        ];
        common::vec_mag(d)
    };
    let v2_mag = common::vec_mag(sol.v2);
    assert!(
        vel_err < tol * v2_mag,
        "{label}: velocity error {vel_err:.2e} (relative {:.2e}), tol={tol:.0e}",
        vel_err / v2_mag
    );
}

// Fixed-geometry multi-rev tests: 90-degree transfer on a near-circular
// orbit with TOF long enough for 1-rev solutions.
// T_circ(r=1, mu=1) = 2*pi ≈ 6.28.  For nrev=1 we need TOF > ~1 period.
// Using TOF ≈ 1.5 * T_circ gives plenty of margin.

#[test]
fn multi_rev_1_long_period_prograde() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "multi_rev_1_long_prograde",
    );
}

#[test]
fn multi_rev_1_short_period_prograde() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "multi_rev_1_short_prograde",
    );
}

#[test]
fn multi_rev_1_long_period_retrograde() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Retrograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "multi_rev_1_long_retrograde",
    );
}

#[test]
fn multi_rev_1_short_period_retrograde() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Retrograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "multi_rev_1_short_retrograde",
    );
}

// 3D geometry with different radii -- nrev=1
#[test]
fn multi_rev_1_3d_long() {
    let t_circ = 2.0 * PI * (2.0_f64.powi(3) / 1.0).sqrt(); // r_ref ≈ 2
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 0.8, 0.6],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "multi_rev_1_3d_long",
    );
}

#[test]
fn multi_rev_1_3d_short() {
    let t_circ = 2.0 * PI * (2.0_f64.powi(3) / 1.0).sqrt();
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 0.8, 0.6],
        1.5 * t_circ,
        1.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "multi_rev_1_3d_short",
    );
}

// Physical units: LEO-to-GEO, nrev=1
#[test]
fn multi_rev_1_leo_geo_long() {
    let mu = 398600.4418_f64;
    let r_ref = (6678.0_f64 + 42164.0) / 2.0;
    let t_circ = 2.0 * PI * (r_ref.powi(3) / mu).sqrt();
    multi_rev_round_trip(
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        2.0 * t_circ,
        mu,
        1,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "multi_rev_1_leo_geo_long",
    );
}

#[test]
fn multi_rev_1_leo_geo_short() {
    let mu = 398600.4418_f64;
    let r_ref = (6678.0_f64 + 42164.0) / 2.0;
    let t_circ = 2.0 * PI * (r_ref.powi(3) / mu).sqrt();
    multi_rev_round_trip(
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        2.0 * t_circ,
        mu,
        1,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "multi_rev_1_leo_geo_short",
    );
}

// nrev=2: higher revolution count.  TOF must exceed ~2 orbital periods.
#[test]
fn multi_rev_2_long_period() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        2.5 * t_circ,
        1.0,
        2,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "multi_rev_2_long",
    );
}

#[test]
fn multi_rev_2_short_period() {
    let t_circ = 2.0 * PI;
    multi_rev_round_trip(
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        2.5 * t_circ,
        1.0,
        2,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "multi_rev_2_short",
    );
}

// -- Randomized multi-rev round-trip tests ------------------------------------
//
// Same pattern as randomized_round_trips but for nrev >= 1.  Uses generous
// TOFs (well above nrev * T_circ) so solutions reliably exist.

/// Run `n` randomized multi-rev round-trip cases.
fn randomized_multi_rev_round_trips(
    seed: &mut u64,
    n: usize,
    mu: f64,
    r_lo: f64,
    r_hi: f64,
    nrev: u32,
    dir: Direction,
    period: MultiRevPeriod,
    tol: f64,
    label: &str,
) -> usize {
    const NEAR_COLLINEAR: f64 = 0.9962;

    let mut tested = 0;
    let mut attempts = 0;

    while tested < n {
        attempts += 1;
        assert!(
            attempts < n * 200,
            "{label}: gave up after {attempts} attempts (only {tested} tested); seed={seed}"
        );

        let r1 = common::rand_r(seed, r_lo, r_hi);
        let r2 = common::rand_r(seed, r_lo, r_hi);

        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        // TOF: nrev+0.5 to nrev+1.5 orbital periods of the larger orbit.
        // This ensures the TOF is above the multi-rev minimum.
        let r1_mag = common::vec_mag(r1);
        let r2_mag = common::vec_mag(r2);
        let r_max = r1_mag.max(r2_mag);
        let t_circ = 2.0 * PI * (r_max * r_max * r_max / mu).sqrt();
        let lo_mult = nrev as f64 + 0.5;
        let hi_mult = nrev as f64 + 1.5;
        let tof = common::rand_f64(seed, lo_mult * t_circ, hi_mult * t_circ);

        let sol = match lambert(mu, r1, r2, tof, nrev, dir, period) {
            Ok(s) => s,
            Err(_) => continue,
        };

        let (r2_prop, v2_prop) = match kepler_propagate(r1, sol.v1, tof, mu) {
            Some(rv) => rv,
            None => continue,
        };

        let pos_err = {
            let d = [r2[0] - r2_prop[0], r2[1] - r2_prop[1], r2[2] - r2_prop[2]];
            common::vec_mag(d)
        };
        assert!(
            pos_err < tol * r2_mag,
            "{label}: position error {pos_err:.2e} (relative {:.2e}), tol={tol:.0e}\n  \
             seed={seed}, attempt={attempts}\n  \
             r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}\n  v1={:.6?}",
            pos_err / r2_mag,
            sol.v1
        );

        let vel_err = {
            let d = [
                sol.v2[0] - v2_prop[0],
                sol.v2[1] - v2_prop[1],
                sol.v2[2] - v2_prop[2],
            ];
            common::vec_mag(d)
        };
        let v2_mag = common::vec_mag(sol.v2);
        assert!(
            vel_err < tol * v2_mag,
            "{label}: velocity error {vel_err:.2e} (relative {:.2e}), tol={tol:.0e}\n  \
             seed={seed}, attempt={attempts}\n  \
             r1={r1:.6?}\n  r2={r2:.6?}\n  tof={tof:.6}\n  v1={:.6?}",
            vel_err / v2_mag,
            sol.v1
        );

        tested += 1;
    }

    tested
}

#[test]
fn random_multi_rev_1_long_50() {
    let mut seed = common::make_seed();
    eprintln!("random_multi_rev_1_long_50: seed = {seed}");

    let tested = randomized_multi_rev_round_trips(
        &mut seed,
        50,
        1.0,
        0.5,
        3.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "random_multi_rev_1_long_50",
    );

    eprintln!("random_multi_rev_1_long_50: {tested} cases passed (seed={seed})");
}

#[test]
fn random_multi_rev_1_short_50() {
    let mut seed = common::make_seed();
    eprintln!("random_multi_rev_1_short_50: seed = {seed}");

    let tested = randomized_multi_rev_round_trips(
        &mut seed,
        50,
        1.0,
        0.5,
        3.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "random_multi_rev_1_short_50",
    );

    eprintln!("random_multi_rev_1_short_50: {tested} cases passed (seed={seed})");
}

#[test]
fn random_multi_rev_1_retrograde_50() {
    let mut seed = common::make_seed();
    eprintln!("random_multi_rev_1_retrograde_50: seed = {seed}");

    // Mix long and short period for retrograde
    let long = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        1.0,
        0.5,
        3.0,
        1,
        Direction::Retrograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "random_multi_rev_1_retro_long",
    );
    let short = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        1.0,
        0.5,
        3.0,
        1,
        Direction::Retrograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "random_multi_rev_1_retro_short",
    );

    eprintln!(
        "random_multi_rev_1_retrograde_50: {long} long + {short} short passed (seed={seed})"
    );
}

#[test]
fn random_multi_rev_2_mixed_50() {
    let mut seed = common::make_seed();
    eprintln!("random_multi_rev_2_mixed_50: seed = {seed}");

    let long = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        1.0,
        0.5,
        3.0,
        2,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "random_multi_rev_2_long",
    );
    let short = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        1.0,
        0.5,
        3.0,
        2,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "random_multi_rev_2_short",
    );

    eprintln!(
        "random_multi_rev_2_mixed_50: {long} long + {short} short passed (seed={seed})"
    );
}

#[test]
fn random_multi_rev_physical_units_50() {
    let mut seed = common::make_seed();
    eprintln!("random_multi_rev_physical_units_50: seed = {seed}");

    let mu = 398600.4418;

    let long = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        mu,
        6400.0,
        30000.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::LongPeriod,
        1e-9,
        "random_multi_rev_phys_long",
    );
    let short = randomized_multi_rev_round_trips(
        &mut seed,
        25,
        mu,
        6400.0,
        30000.0,
        1,
        Direction::Prograde,
        MultiRevPeriod::ShortPeriod,
        1e-9,
        "random_multi_rev_phys_short",
    );

    eprintln!(
        "random_multi_rev_physical_units_50: {long} long + {short} short passed (seed={seed})"
    );
}
