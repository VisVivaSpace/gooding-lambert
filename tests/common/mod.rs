//! Shared test utilities: seeded PRNG and geometry helpers.
//!
//! Included via `mod common;` in each integration test that needs random geometry.
#![allow(dead_code)]

use std::f64::consts::PI;

// ── Seeded PRNG (xorshift64) ─────────────────────────────────────────────────

pub fn make_seed() -> u64 {
    let d = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap();
    let s = d
        .as_secs()
        .wrapping_mul(6364136223846793005)
        .wrapping_add(d.subsec_nanos() as u64);
    if s == 0 { 1 } else { s }
}

pub fn xorshift(s: &mut u64) -> u64 {
    *s ^= *s << 13;
    *s ^= *s >> 7;
    *s ^= *s << 17;
    *s
}

pub fn rand_f64(s: &mut u64, lo: f64, hi: f64) -> f64 {
    let bits = xorshift(s) >> 11; // 53 random bits
    lo + (bits as f64 / (1u64 << 53) as f64) * (hi - lo)
}

/// Random 3D position vector: random direction, magnitude in [r_lo, r_hi].
pub fn rand_r(s: &mut u64, r_lo: f64, r_hi: f64) -> [f64; 3] {
    let r = rand_f64(s, r_lo, r_hi);
    let theta = rand_f64(s, 0.0, PI);
    let phi = rand_f64(s, 0.0, 2.0 * PI);
    [
        r * theta.sin() * phi.cos(),
        r * theta.sin() * phi.sin(),
        r * theta.cos(),
    ]
}

// ── Geometry helpers ─────────────────────────────────────────────────────────

pub fn vec_mag(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

pub fn cos_transfer_angle(r1: [f64; 3], r2: [f64; 3]) -> f64 {
    (r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) / (vec_mag(r1) * vec_mag(r2))
}
