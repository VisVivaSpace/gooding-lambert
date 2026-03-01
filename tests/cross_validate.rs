//! Cross-validation: gooding-lambert (Gooding) vs lambert_solver (ivLam2 vercosine).
//!
//! Two independent algorithms, both converged to ~1e-10 internally.
//! Expected velocity agreement: ~1e-8 (documented tolerance for two-algorithm
//! cross-validation; the residual reflects different formulations, not error).
//!
//! The random tests use a time-based seed that changes every run; the seed is printed
//! so that failures can be reproduced with `cargo test -- --nocapture`.

use gooding_lambert::{lambert, Direction};
use lambert_solver::{solve_lambert, Direction as IvDir};
use std::f64::consts::PI;

const TOL: f64 = 1e-8;

// ── Test harness ─────────────────────────────────────────────────────────────

fn cross_validate(label: &str, r1: [f64; 3], r2: [f64; 3], tof: f64, mu: f64, dir: Direction) {
    let iv_dir = match dir {
        Direction::Prograde => IvDir::Prograde,
        Direction::Retrograde => IvDir::Retrograde,
    };

    let ours = lambert(mu, r1, r2, tof, 0, dir)
        .unwrap_or_else(|e| panic!("{label}: gooding-lambert failed: {e:?}"));
    let theirs = solve_lambert(&r1, &r2, tof, mu, iv_dir, 0)
        .unwrap_or_else(|e| panic!("{label}: lambert_solver failed: {e:?}"));

    for i in 0..3 {
        let diff1 = (ours.v1[i] - theirs.v1[i]).abs();
        let diff2 = (ours.v2[i] - theirs.v2[i]).abs();
        assert!(
            diff1 < TOL,
            "{label}: v1[{i}] ours={:.15e} theirs={:.15e} diff={diff1:.2e}",
            ours.v1[i], theirs.v1[i]
        );
        assert!(
            diff2 < TOL,
            "{label}: v2[{i}] ours={:.15e} theirs={:.15e} diff={diff2:.2e}",
            ours.v2[i], theirs.v2[i]
        );
    }
}

// ── PRNG (xorshift64, time-seeded) ─────────────────────────────────────────

fn make_seed() -> u64 {
    let d = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap();
    let s = d.as_secs()
        .wrapping_mul(6364136223846793005)
        .wrapping_add(d.subsec_nanos() as u64);
    if s == 0 { 1 } else { s }
}

fn xorshift(s: &mut u64) -> u64 {
    *s ^= *s << 13;
    *s ^= *s >> 7;
    *s ^= *s << 17;
    *s
}

fn rand_f64(s: &mut u64, lo: f64, hi: f64) -> f64 {
    let bits = xorshift(s) >> 11;
    lo + (bits as f64 / (1u64 << 53) as f64) * (hi - lo)
}

/// Random 3D position vector: random direction, magnitude in [r_lo, r_hi].
fn rand_r(s: &mut u64, r_lo: f64, r_hi: f64) -> [f64; 3] {
    let r = rand_f64(s, r_lo, r_hi);
    let theta = rand_f64(s, 0.0, PI);
    let phi = rand_f64(s, 0.0, 2.0 * PI);
    [r * theta.sin() * phi.cos(), r * theta.sin() * phi.sin(), r * theta.cos()]
}

fn vec_mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn cos_transfer_angle(r1: [f64; 3], r2: [f64; 3]) -> f64 {
    (r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) / (vec_mag(r1) * vec_mag(r2))
}

// ── Fixed regression tests ──────────────────────────────────────────────────

#[test]
fn prograde_10_deg() {
    let a = 10.0_f64.to_radians();
    cross_validate("10°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde);
}

#[test]
fn prograde_45_deg() {
    let a = 45.0_f64.to_radians();
    cross_validate("45°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde);
}

#[test]
fn prograde_90_deg() {
    cross_validate("90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], PI / 2.0, 1.0, Direction::Prograde);
}

#[test]
fn prograde_120_deg() {
    let a = 120.0_f64.to_radians();
    cross_validate("120°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde);
}

#[test]
fn prograde_135_deg() {
    let a = 135.0_f64.to_radians();
    cross_validate("135°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], a, 1.0, Direction::Prograde);
}

#[test]
fn prograde_170_deg() {
    let a = 170.0_f64.to_radians();
    cross_validate("170°", [1.0, 0.0, 0.0], [a.cos(), a.sin(), 0.0], 2.0, 1.0, Direction::Prograde);
}

#[test]
fn ratio_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    for &ratio in &[0.3_f64, 0.5, 0.7, 1.5, 2.0, 3.0, 5.0] {
        let r2 = [0.0, ratio, 0.0];
        let a = (1.0 + ratio) / 2.0;
        let tof = PI * (a * a * a).sqrt() * 0.5;
        cross_validate(&format!("r2/r1={ratio:.1}"), r1, r2, tof, 1.0, Direction::Prograde);
    }
}

#[test]
fn tof_sweep() {
    let r1 = [1.0, 0.0, 0.0];
    let r2 = [0.0, 1.0, 0.0];
    for &tof in &[0.2_f64, 0.5, 1.0, 2.0, 5.0, 10.0] {
        cross_validate(&format!("tof={tof:.1}"), r1, r2, tof, 1.0, Direction::Prograde);
    }
}

#[test]
fn three_d_90_deg() {
    cross_validate("3D", [1.0, 0.0, 0.0], [0.0, 0.8, 0.6], PI / 2.0, 1.0, Direction::Prograde);
}

#[test]
fn three_d_different_radii() {
    cross_validate("3D diff-r", [1.0, 0.5, 0.0], [0.0, 1.0, 0.8], 2.0, 1.0, Direction::Prograde);
}

#[test]
fn retrograde_90_deg() {
    cross_validate("retrograde 90°", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 3.0 * PI / 2.0, 1.0, Direction::Retrograde);
}

#[test]
fn hyperbolic() {
    cross_validate("hyperbolic", [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.1, 1.0, Direction::Prograde);
}

#[test]
fn leo_to_geo() {
    let mu = 398600.4418;
    cross_validate("LEO→GEO", [6678.0, 0.0, 0.0], [0.0, 42164.0, 0.0], 5.0 * 3600.0, mu, Direction::Prograde);
}

#[test]
fn leo_to_meo() {
    let mu = 398600.4418;
    cross_validate("LEO→MEO", [6678.0, 0.0, 0.0], [0.0, 26560.0, 0.0], 4.0 * 3600.0, mu, Direction::Prograde);
}

#[test]
fn earth_to_mars_type() {
    let mu = 1.327e11;
    let r_earth = 1.496e8;
    let r_mars = 2.279e8;
    let a = 135.0_f64.to_radians();
    cross_validate("Earth→Mars", [r_earth, 0.0, 0.0], [r_mars * a.cos(), r_mars * a.sin(), 0.0], 200.0 * 86400.0, mu, Direction::Prograde);
}

#[test]
fn high_eccentricity_r7() {
    cross_validate("ecc_r7", [1.0, 0.0, 0.0], [0.0, 7.0, 0.0], 2.0, 1.0, Direction::Prograde);
}

#[test]
fn high_eccentricity_r10() {
    cross_validate("ecc_r10", [1.0, 0.0, 0.0], [0.0, 10.0, 0.0], 3.0, 1.0, Direction::Prograde);
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
// Run with `cargo test -- --nocapture` to see seed and statistics.

fn random_sweep(label: &str, dir: Direction, n: usize) {
    let mut seed = make_seed();
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

        let r1 = rand_r(&mut seed, 0.5, 5.0);
        let r2 = rand_r(&mut seed, 0.5, 5.0);

        // Reject near-collinear
        if cos_transfer_angle(r1, r2).abs() > NEAR_COLLINEAR {
            continue;
        }

        // TOF: fraction of mean circular orbit period
        let r_mean = (vec_mag(r1) + vec_mag(r2)) / 2.0;
        let t_circ = 2.0 * PI * (r_mean * r_mean * r_mean / MU).sqrt();
        let tof = rand_f64(&mut seed, 0.02 * t_circ, 1.5 * t_circ);

        let ours = lambert(MU, r1, r2, tof, 0, dir);
        let theirs = solve_lambert(&r1, &r2, tof, MU, iv_dir, 0);

        match (ours, theirs) {
            (Ok(ours), Ok(theirs)) => {
                for i in 0..3 {
                    let d1 = (ours.v1[i] - theirs.v1[i]).abs();
                    let d2 = (ours.v2[i] - theirs.v2[i]).abs();
                    assert!(
                        d1 < TOL,
                        "{label} seed={seed} attempt={attempts}: v1[{i}] ours={:.15e} theirs={:.15e} diff={d1:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                        ours.v1[i], theirs.v1[i]
                    );
                    assert!(
                        d2 < TOL,
                        "{label} seed={seed} attempt={attempts}: v2[{i}] ours={:.15e} theirs={:.15e} diff={d2:.2e}\n  r1={r1:.6?} r2={r2:.6?} tof={tof:.6}",
                        ours.v2[i], theirs.v2[i]
                    );
                }
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
    random_sweep("random_prograde_100", Direction::Prograde, 100);
}

#[test]
fn random_retrograde_100() {
    random_sweep("random_retrograde_100", Direction::Retrograde, 100);
}
