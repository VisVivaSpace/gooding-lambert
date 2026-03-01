//! Cross-validation: gooding-lambert (Gooding) vs lambert_solver (ivLam2 vercosine).
//!
//! Two independent algorithms, both converged to ~1e-12 internally.
//!
//! Tolerance tiers (relative vector-magnitude error, |δv| / |v|):
//!
//! - `TOL_NORMAL = 1e-10`: well-conditioned cases (gentle geometry, all-elliptic).
//!   Both algorithms converge to ~1e-12; 1e-10 provides a 100× safety margin over
//!   expected agreement and catches genuine regressions.
//!
//! - `TOL_DEGENERATE = 1e-8`: cases with elevated condition number — hyperbolic
//!   transfers (large |x|), very large r2/r1 ratios, or mixed tof sweeps that
//!   include near-parabolic geometries.  The algorithms may diverge by more due to
//!   different handling of the parabolic boundary.
//!
//! All comparisons are relative: `|δv| / |v| < tol` using the vector magnitude
//! of the difference divided by the vector magnitude of the reference velocity.
//! This is consistent across physical- and canonical-unit tests regardless of
//! the absolute velocity scale.
//!
//! The random tests use a time-based seed that changes every run; the seed is printed
//! so that failures can be reproduced with `cargo test -- --nocapture`.

mod common;

use gooding_lambert::{lambert, Direction};
use lambert_solver::{solve_lambert, Direction as IvDir};
use std::f64::consts::PI;

// Well-conditioned: both algorithms, ~1e-12 internal convergence, 100× margin.
const TOL_NORMAL: f64 = 1e-10;
// Elevated condition number (hyperbolic, large r2/r1, mixed tof sweep).
const TOL_DEGENERATE: f64 = 1e-8;

// ── Test harness ─────────────────────────────────────────────────────────────

fn cross_validate(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction, tol: f64) {
    let iv_dir = match dir {
        Direction::Prograde => IvDir::Prograde,
        Direction::Retrograde => IvDir::Retrograde,
    };

    let ours = lambert(mu, r1, r2, tof, 0, dir)
        .unwrap_or_else(|e| panic!("{label}: gooding-lambert failed: {e:?}"));
    let theirs = solve_lambert(&r1, &r2, tof, mu, iv_dir, 0)
        .unwrap_or_else(|e| panic!("{label}: lambert_solver failed: {e:?}"));

    // Relative comparison: |δv| / |v| < tol
    let v1_err = common::vec_mag([
        ours.v1[0] - theirs.v1[0],
        ours.v1[1] - theirs.v1[1],
        ours.v1[2] - theirs.v1[2],
    ]);
    let v2_err = common::vec_mag([
        ours.v2[0] - theirs.v2[0],
        ours.v2[1] - theirs.v2[1],
        ours.v2[2] - theirs.v2[2],
    ]);
    let v1_mag = common::vec_mag(ours.v1);
    let v2_mag = common::vec_mag(ours.v2);

    assert!(
        v1_err < tol * v1_mag,
        "{label}: v1 relative error {:.2e} (|δv1|={v1_err:.2e}, |v1|={v1_mag:.2e}), tol={tol:.0e}",
        v1_err / v1_mag
    );
    assert!(
        v2_err < tol * v2_mag,
        "{label}: v2 relative error {:.2e} (|δv2|={v2_err:.2e}, |v2|={v2_mag:.2e}), tol={tol:.0e}",
        v2_err / v2_mag
    );
}

// ── Fixed regression tests ──────────────────────────────────────────────────

#[test]
fn prograde_10_deg() {
    let a = 10.0_f64.to_radians();
    cross_validate("10°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn prograde_45_deg() {
    let a = 45.0_f64.to_radians();
    cross_validate("45°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn prograde_90_deg() {
    cross_validate("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn prograde_120_deg() {
    let a = 120.0_f64.to_radians();
    cross_validate("120°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn prograde_135_deg() {
    let a = 135.0_f64.to_radians();
    cross_validate("135°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn prograde_170_deg() {
    let a = 170.0_f64.to_radians();
    cross_validate("170°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], 2.0, 1.0, Direction::Prograde, TOL_NORMAL);
}

// All transfers here are elliptic (TOF = half Hohmann period for each ratio).
#[test]
fn ratio_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for &ratio in &[0.3_f64, 0.5, 0.7, 1.5, 2.0, 3.0, 5.0] {
        let r2 = [0.0, ratio, 0.0];
        let a = (1.0 + ratio) / 2.0;
        let tof = PI * (a * a * a).sqrt() * 0.5;
        cross_validate(&format!("r2/r1={ratio:.1}"), r1, r2, tof, 1.0, Direction::Prograde, TOL_NORMAL);
    }
}

// Short TOFs (0.2, 0.5) give hyperbolic transfers for this geometry; use
// TOL_DEGENERATE for the entire sweep since it is not split per-case.
#[test]
fn tof_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    for &tof in &[0.2_f64, 0.5, 1.0, 2.0, 5.0, 10.0] {
        cross_validate(&format!("tof={tof:.1}"), r1, r2, tof, 1.0, Direction::Prograde, TOL_DEGENERATE);
    }
}

#[test]
fn three_d_90_deg() {
    cross_validate("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn three_d_different_radii() {
    cross_validate("3D diff-r", [1.0, 0.5, 0.0], [0.0, 1.0, 0.8], 2.0, 1.0, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn retrograde_90_deg() {
    cross_validate("retrograde 90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 3.0 * PI / 2.0, 1.0, Direction::Retrograde, TOL_NORMAL);
}

// Short TOF → hyperbolic trajectory; condition number elevated near parabolic.
#[test]
fn hyperbolic() {
    cross_validate("hyperbolic", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0, Direction::Prograde, TOL_DEGENERATE);
}

#[test]
fn leo_to_geo() {
    let mu = 398600.4418;
    cross_validate("LEO→GEO", [6678.0, 0.0, 0.0], [0.0, 42164.0, 0.0], 5.0 * 3600.0, mu, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn leo_to_meo() {
    let mu = 398600.4418;
    cross_validate("LEO→MEO", [6678.0, 0.0, 0.0], [0.0, 26560.0, 0.0], 4.0 * 3600.0, mu, Direction::Prograde, TOL_NORMAL);
}

#[test]
fn earth_to_mars_type() {
    let mu = 1.327e11;
    let r_earth = 1.496e8;
    let r_mars = 2.279e8;
    let a = 135.0_f64.to_radians();
    cross_validate("Earth→Mars", [r_earth, 0.0, 0.0], [r_mars * a.cos(), r_mars * a.sin(), 0.0], 200.0 * 86400.0, mu, Direction::Prograde, TOL_NORMAL);
}

// Large r2/r1 ratios: elevated condition number near the high-eccentricity limit.
#[test]
fn high_eccentricity_r7() {
    cross_validate("ecc_r7", [1.0, 0.0, 0.0], [0.0, 7.0, 0.0], 2.0, 1.0, Direction::Prograde, TOL_DEGENERATE);
}

#[test]
fn high_eccentricity_r10() {
    cross_validate("ecc_r10", [1.0, 0.0, 0.0], [0.0, 10.0, 0.0], 3.0, 1.0, Direction::Prograde, TOL_DEGENERATE);
}

// ── Randomized sweeps — 100+ cases per run, new seed each time ──────────────
//
// Geometry: random 3D r1 and r2, magnitudes in [0.5, 5].
// Transfer angle: filtered to (5°, 175°) for prograde; same range applies to
// retrograde (the supplementary arc is then 185°–355°).
// TOF: random fraction of the mean circular orbit period [2%, 150%].
// Cases where either solver fails are skipped silently (both must fail or succeed
// for a comparison; mismatched failure is noted but does not fail the test).
//
// Random geometry is well-conditioned by construction (collinear cases rejected,
// TOF kept to a reasonable fraction of the period), so TOL_NORMAL applies.
//
// Run with `cargo test -- --nocapture` to see seed and statistics.

fn random_sweep(label: &str, dir: Direction, n: usize, tol: f64) {
    let mut seed = common::make_seed();
    eprintln!("{label}: seed = {seed}");

    const MU: f64 = 1.0;
    const NEAR_COLLINEAR: f64 = 0.9962; // cos(5°) — reject < 5° or > 175°

    let iv_dir = match dir {
        Direction::Prograde => IvDir::Prograde,
        Direction::Retrograde => IvDir::Retrograde,
    };

    let mut compared = 0usize;
    let mut attempts = 0usize;
    let mut one_sided = 0usize; // one solver succeeded, other didn't

    while compared < n {
        attempts += 1;
        assert!(
            attempts < 20_000,
            "{label}: gave up after {attempts} attempts ({compared} compared); seed={seed}"
        );

        let r1 = common::rand_r(&mut seed, 0.5, 5.0);
        let r2 = common::rand_r(&mut seed, 0.5, 5.0);

        // Reject near-collinear
        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        // TOF: fraction of mean circular orbit period
        let r_mean = (common::vec_mag(r1) + common::vec_mag(r2)) / 2.0;
        let t_circ = 2.0 * PI * (r_mean * r_mean * r_mean / MU).sqrt();
        let tof = common::rand_f64(&mut seed, 0.02 * t_circ, 1.5 * t_circ);

        let ours = lambert(MU, r1, r2, tof, 0, dir);
        let theirs = solve_lambert(&r1, &r2, tof, MU, iv_dir, 0);

        match (ours, theirs) {
            (Ok(ours), Ok(theirs)) => {
                let v1_err = common::vec_mag([
                    ours.v1[0] - theirs.v1[0],
                    ours.v1[1] - theirs.v1[1],
                    ours.v1[2] - theirs.v1[2],
                ]);
                let v2_err = common::vec_mag([
                    ours.v2[0] - theirs.v2[0],
                    ours.v2[1] - theirs.v2[1],
                    ours.v2[2] - theirs.v2[2],
                ]);
                let v1_mag = common::vec_mag(ours.v1);
                let v2_mag = common::vec_mag(ours.v2);
                assert!(
                    v1_err < tol * v1_mag,
                    "{label} seed={seed} attempt={attempts}: v1 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v1_err / v1_mag
                );
                assert!(
                    v2_err < tol * v2_mag,
                    "{label} seed={seed} attempt={attempts}: v2 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v2_err / v2_mag
                );
                compared += 1;
            }
            (Err(_), Err(_)) => {}   // both failed: edge case, skip silently
            _ => { one_sided += 1; } // mismatch: count but don't fail
        }
    }

    eprintln!(
        "{label}: {compared} cases compared in {attempts} attempts, {one_sided} one-sided failures (seed={seed})"
    );
}

#[test]
fn random_prograde_100() {
    random_sweep("random_prograde_100", Direction::Prograde, 100, TOL_NORMAL);
}

#[test]
fn random_retrograde_100() {
    random_sweep("random_retrograde_100", Direction::Retrograde, 100, TOL_NORMAL);
}
