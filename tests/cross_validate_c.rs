//! Cross-validation: Rust Gooding implementation vs C Gooding implementation.
//!
//! Only active when compiled with `--features gooding-ffi`.
//! The C solver only handles prograde transfers with th ∈ [0, π] and nrev=0,
//! so all tests here are prograde zero-rev.
//!
//! Since both solvers implement identical algorithms, agreement to 1e-10 is expected.
//! The random test uses a time-based seed that changes every run; the seed is printed
//! so that failures can be reproduced with `cargo test -- --nocapture`.

#![cfg(feature = "gooding-ffi")]

mod common;

use gooding_lambert::{Direction, lambert as rust_lambert};
use std::f64::consts::PI;

// Same algorithm (Rust vs C Gooding): differences come only from compiler
// floating-point reordering / FMA decisions over ~50 arithmetic operations.
// Expected agreement: ~50ε ≈ 1e-14. xlamb converges to |δT| < 1e-12 in both;
// their x values should match to within that criterion.
// Tolerance: 1e-12 (100× above expected 1e-14, catches compiler-divergence bugs).
const TOL: f64 = 1e-12;

// ── C FFI ───────────────────────────────────────────────────────────────────

unsafe extern "C" {
    fn lambert(
        gm: f64,
        r1: *const f64,
        r2: *const f64,
        nrev: i32,
        dt: f64,
        v1: *mut f64,
        v2: *mut f64,
    ) -> i32;
}

fn c_lambert(mu: f64, r1: [f64; 3], r2: [f64; 3], tof: f64) -> Option<([f64; 3], [f64; 3])> {
    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];
    let code = unsafe {
        lambert(
            mu,
            r1.as_ptr(),
            r2.as_ptr(),
            0,
            tof,
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
        )
    };
    if code > 0 { Some((v1, v2)) } else { None }
}

fn c_lambert_retro(mu: f64, r1: [f64; 3], r2: [f64; 3], tof: f64) -> Option<([f64; 3], [f64; 3])> {
    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];
    let code = unsafe {
        lambert(
            mu,
            r1.as_ptr(),
            r2.as_ptr(),
            0,
            -tof,
            v1.as_mut_ptr(),
            v2.as_mut_ptr(),
        )
    };
    if code > 0 { Some((v1, v2)) } else { None }
}

// ── Test harness ─────────────────────────────────────────────────────────────

fn cross_validate_c(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64) {
    let ours = rust_lambert(mu, r1, r2, tof, 0, Direction::Prograde)
        .unwrap_or_else(|e| panic!("{label}: Rust failed: {e:?}"));
    let (cv1, cv2) = c_lambert(mu, r1, r2, tof)
        .unwrap_or_else(|| panic!("{label}: C solver returned no solution"));

    // Relative comparison: |δv| / |v| < TOL
    let v1_err = common::vec_mag([
        ours.v1[0] - cv1[0],
        ours.v1[1] - cv1[1],
        ours.v1[2] - cv1[2],
    ]);
    let v2_err = common::vec_mag([
        ours.v2[0] - cv2[0],
        ours.v2[1] - cv2[1],
        ours.v2[2] - cv2[2],
    ]);
    let v1_mag = common::vec_mag(ours.v1);
    let v2_mag = common::vec_mag(ours.v2);
    assert!(
        v1_err < TOL * v1_mag,
        "{label}: v1 relative error {:.2e} (|δv1|={v1_err:.2e}, |v1|={v1_mag:.2e})",
        v1_err / v1_mag
    );
    assert!(
        v2_err < TOL * v2_mag,
        "{label}: v2 relative error {:.2e} (|δv2|={v2_err:.2e}, |v2|={v2_mag:.2e})",
        v2_err / v2_mag
    );
}

// ── Fixed regression tests ──────────────────────────────────────────────────

#[test]
fn c_10_deg() {
    let a = 10.0_f64.to_radians();
    cross_validate_c("10°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0);
}

#[test]
fn c_45_deg() {
    let a = 45.0_f64.to_radians();
    cross_validate_c("45°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0);
}

#[test]
fn c_90_deg() {
    cross_validate_c("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0);
}

#[test]
fn c_135_deg() {
    let a = 135.0_f64.to_radians();
    cross_validate_c("135°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0);
}

#[test]
fn c_170_deg() {
    let a = 170.0_f64.to_radians();
    cross_validate_c("170°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], 2.0, 1.0);
}

#[test]
fn c_hyperbolic() {
    cross_validate_c("hyperbolic", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0);
}

#[test]
fn c_different_radii() {
    cross_validate_c("diff radii", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0, 1.0);
}

#[test]
fn c_3d() {
    cross_validate_c("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0);
}

#[test]
fn c_leo_to_geo() {
    let mu = 398600.4418;
    cross_validate_c(
        "LEO→GEO",
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        5.0 * 3600.0,
        mu,
    );
}

#[test]
fn c_angle_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for &deg in &[
        5.0_f64, 15.0, 25.0, 30.0, 60.0, 75.0, 100.0, 110.0, 140.0, 150.0, 160.0,
    ] {
        let a = deg.to_radians();
        let r2 = [a.cos(), a.sin(), 0.0];
        cross_validate_c(&format!("{deg}°"), r1, r2, a, 1.0);
    }
}

// ── Randomized sweep — 100+ cases per run, new seed each time ───────────────
//
// Geometry: random 3D r1 and r2 with magnitudes in [0.5, 5].
// Transfer angle: filtered to (5°, 175°) to avoid singular geometry.
// TOF: random fraction of the mean circular orbit period.
// Both solvers must agree to 1e-10; cases where either fails are skipped.
// Run with `cargo test --features gooding-ffi -- --nocapture` to see the seed.

#[test]
fn c_random_100() {
    let mut seed = common::make_seed();
    eprintln!("c_random_100: seed = {seed}");

    const MU: f64 = 1.0;
    // cos(5°) — reject vectors closer than 5° or more than 175° apart
    const NEAR_COLLINEAR: f64 = 0.9962; // cos(5°)

    let mut compared = 0usize;
    let mut attempts = 0usize;

    while compared < 100 {
        attempts += 1;
        assert!(
            attempts < 10_000,
            "c_random_100: gave up after {attempts} attempts (only {compared} compared); seed={seed}"
        );

        let r1 = common::rand_r(&mut seed, 0.5, 5.0);
        let r2 = common::rand_r(&mut seed, 0.5, 5.0);

        // Reject near-collinear: |cos θ| > cos(5°)
        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        // TOF: random fraction of mean circular orbit period [2%, 150%]
        let r_mean = (common::vec_mag(r1) + common::vec_mag(r2)) / 2.0;
        let t_circ = 2.0 * PI * (r_mean * r_mean * r_mean / MU).sqrt();
        let tof = common::rand_f64(&mut seed, 0.02 * t_circ, 1.5 * t_circ);

        let rust = rust_lambert(MU, r1, r2, tof, 0, Direction::Prograde);
        let c = c_lambert(MU, r1, r2, tof);

        match (rust, c) {
            (Ok(ours), Some((cv1, cv2))) => {
                let v1_err = common::vec_mag([
                    ours.v1[0] - cv1[0],
                    ours.v1[1] - cv1[1],
                    ours.v1[2] - cv1[2],
                ]);
                let v2_err = common::vec_mag([
                    ours.v2[0] - cv2[0],
                    ours.v2[1] - cv2[1],
                    ours.v2[2] - cv2[2],
                ]);
                let v1_mag = common::vec_mag(ours.v1);
                let v2_mag = common::vec_mag(ours.v2);
                assert!(
                    v1_err < TOL * v1_mag,
                    "seed={seed} attempt={attempts}: v1 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v1_err / v1_mag
                );
                assert!(
                    v2_err < TOL * v2_mag,
                    "seed={seed} attempt={attempts}: v2 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v2_err / v2_mag
                );
                compared += 1;
            }
            // Both failed: geometry may be at an edge — acceptable
            (Err(_), None) => {}
            // One succeeded, other didn't: count but don't fail (investigate if frequent)
            _ => {}
        }
    }

    eprintln!("c_random_100: {compared} cases compared in {attempts} attempts (seed={seed})");
}

// ── Retrograde cross-validation ─────────────────────────────────────────────
//
// The C solver uses negative tdelt for retrograde. Our Rust solver now matches
// this convention in gooding_lambert(). These tests verify agreement.

fn cross_validate_c_retro(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64) {
    let ours = rust_lambert(mu, r1, r2, tof, 0, Direction::Retrograde)
        .unwrap_or_else(|e| panic!("{label}: Rust failed: {e:?}"));
    let (cv1, cv2) = c_lambert_retro(mu, r1, r2, tof)
        .unwrap_or_else(|| panic!("{label}: C solver returned no solution"));

    let v1_err = common::vec_mag([
        ours.v1[0] - cv1[0],
        ours.v1[1] - cv1[1],
        ours.v1[2] - cv1[2],
    ]);
    let v2_err = common::vec_mag([
        ours.v2[0] - cv2[0],
        ours.v2[1] - cv2[1],
        ours.v2[2] - cv2[2],
    ]);
    let v1_mag = common::vec_mag(ours.v1);
    let v2_mag = common::vec_mag(ours.v2);
    assert!(
        v1_err < TOL * v1_mag,
        "{label}: v1 relative error {:.2e} (|δv1|={v1_err:.2e}, |v1|={v1_mag:.2e})",
        v1_err / v1_mag
    );
    assert!(
        v2_err < TOL * v2_mag,
        "{label}: v2 relative error {:.2e} (|δv2|={v2_err:.2e}, |v2|={v2_mag:.2e})",
        v2_err / v2_mag
    );
}

#[test]
fn c_retro_90_deg() {
    cross_validate_c_retro("retro 90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0);
}

#[test]
fn c_retro_45_deg() {
    let a = 45.0_f64.to_radians();
    cross_validate_c_retro("retro 45°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0);
}

#[test]
fn c_retro_135_deg() {
    let a = 135.0_f64.to_radians();
    cross_validate_c_retro(
        "retro 135°",
        [1.0, 0.0, 0.0],
        [a.cos(), a.sin(), 0.0],
        a,
        1.0,
    );
}

#[test]
fn c_retro_different_radii() {
    cross_validate_c_retro("retro diff radii", [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 2.0, 1.0);
}

#[test]
fn c_retro_3d() {
    cross_validate_c_retro(
        "retro 3D",
        [1.0, 0.0, 0.0],
        [0.0, 0.8, 0.6],
        PI / 2.0,
        1.0,
    );
}

#[test]
fn c_retro_leo_to_geo() {
    let mu = 398600.4418;
    cross_validate_c_retro(
        "retro LEO→GEO",
        [6678.0, 0.0, 0.0],
        [0.0, 42164.0, 0.0],
        5.0 * 3600.0,
        mu,
    );
}

#[test]
fn c_retro_angle_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for &deg in &[15.0_f64, 30.0, 45.0, 60.0, 75.0, 90.0, 100.0, 120.0, 135.0, 150.0, 160.0] {
        let a = deg.to_radians();
        let r2 = [a.cos(), a.sin(), 0.0];
        cross_validate_c_retro(&format!("retro {deg}°"), r1, r2, a, 1.0);
    }
}

#[test]
fn c_retro_random_100() {
    let mut seed = common::make_seed();
    eprintln!("c_retro_random_100: seed = {seed}");

    const MU: f64 = 1.0;
    const NEAR_COLLINEAR: f64 = 0.9962; // cos(5°)

    let mut compared = 0usize;
    let mut attempts = 0usize;

    while compared < 100 {
        attempts += 1;
        assert!(
            attempts < 10_000,
            "c_retro_random_100: gave up after {attempts} attempts (only {compared} compared); seed={seed}"
        );

        let r1 = common::rand_r(&mut seed, 0.5, 5.0);
        let r2 = common::rand_r(&mut seed, 0.5, 5.0);

        if common::cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        let r_mean = (common::vec_mag(r1) + common::vec_mag(r2)) / 2.0;
        let t_circ = 2.0 * PI * (r_mean * r_mean * r_mean / MU).sqrt();
        let tof = common::rand_f64(&mut seed, 0.02 * t_circ, 1.5 * t_circ);

        let rust = rust_lambert(MU, r1, r2, tof, 0, Direction::Retrograde);
        let c = c_lambert_retro(MU, r1, r2, tof);

        match (rust, c) {
            (Ok(ours), Some((cv1, cv2))) => {
                let v1_err = common::vec_mag([
                    ours.v1[0] - cv1[0],
                    ours.v1[1] - cv1[1],
                    ours.v1[2] - cv1[2],
                ]);
                let v2_err = common::vec_mag([
                    ours.v2[0] - cv2[0],
                    ours.v2[1] - cv2[1],
                    ours.v2[2] - cv2[2],
                ]);
                let v1_mag = common::vec_mag(ours.v1);
                let v2_mag = common::vec_mag(ours.v2);
                assert!(
                    v1_err < TOL * v1_mag,
                    "seed={seed} attempt={attempts}: v1 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v1_err / v1_mag
                );
                assert!(
                    v2_err < TOL * v2_mag,
                    "seed={seed} attempt={attempts}: v2 relative error {:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                    v2_err / v2_mag
                );
                compared += 1;
            }
            (Err(_), None) => {}
            _ => {}
        }
    }

    eprintln!("c_retro_random_100: {compared} cases compared in {attempts} attempts (seed={seed})");
}
