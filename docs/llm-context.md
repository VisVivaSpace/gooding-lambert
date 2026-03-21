# LLM Context for gooding-lambert

## What This Crate Does

Solves Lambert's problem: given two position vectors and a time of flight, compute the initial
and final velocity vectors of the connecting Keplerian orbit. Uses Gooding's (1990) method, which
is robust across all orbit types (elliptic, parabolic, hyperbolic), handles multi-revolution
solutions, and avoids the singularities common in other formulations.

## Key Types and Functions

```rust
use gooding_lambert::{lambert, Direction, LambertSolution, LambertError};

// Primary entry point
fn lambert(
    mu: f64,          // gravitational parameter (same units as r³/t²)
    r1: [f64; 3],     // departure position vector
    r2: [f64; 3],     // arrival position vector
    tof: f64,         // time of flight (positive)
    nrev: u32,        // complete revolutions before arrival (0 = single arc)
    dir: Direction,   // Prograde or Retrograde
) -> Result<LambertSolution, LambertError>

// Solution
struct LambertSolution {
    v1: [f64; 3],  // departure velocity vector
    v2: [f64; 3],  // arrival velocity vector
}

// Transfer direction
enum Direction {
    Prograde,    // transfer angle ∈ (0°, 180°); angular momentum parallel to r1 × r2
    Retrograde,  // transfer angle ∈ (180°, 360°); angular momentum parallel to r2 × r1
}

// Errors
enum LambertError {
    SingularTransfer,         // 180° transfer — plane undefined
    NoSolution,               // TOF below minimum for requested nrev
    ConvergenceFailed,        // Householder iteration did not converge
    InvalidInput(&'static str), // zero radius, non-positive TOF, NaN, etc.
}
```

## Common Usage Patterns

### Prograde single-revolution transfer (most common)

```rust
use gooding_lambert::{lambert, Direction};

let mu = 398600.4418_f64;             // Earth GM (km³/s²)
let r1 = [6678.0, 0.0, 0.0];         // LEO departure (km)
let r2 = [0.0, 42164.0, 0.0];        // GEO arrival (km)
let tof = 5.0 * 3600.0;              // 5 hours (seconds)

let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde).unwrap();
// sol.v1 = departure velocity vector (km/s)
// sol.v2 = arrival velocity vector (km/s)
```

### Retrograde transfer

```rust
// Retrograde uses the long-arc (> 180°); tof must be larger accordingly
let sol = lambert(mu, r1, r2, tof, 0, Direction::Retrograde)?;
```

### Multi-revolution transfer

```rust
// nrev=1: one full loop before arrival; tof must exceed the minimum for 1 rev
// Returns the long-period solution for nrev > 0
let sol = lambert(mu, r1, r2, tof, 1, Direction::Prograde)?;
```

### Handling errors

```rust
match lambert(mu, r1, r2, tof, nrev, dir) {
    Ok(sol) => { /* use sol.v1, sol.v2 */ }
    Err(LambertError::NoSolution) => { /* tof too short for this nrev */ }
    Err(LambertError::SingularTransfer) => { /* r1 and r2 are anti-parallel */ }
    Err(e) => { /* other error */ }
}
```

## Unit Conventions

All inputs must be in consistent units — the crate does not convert units.

| Unit system | mu        | r      | tof     | v output  |
|-------------|-----------|--------|---------|-----------|
| SI          | m³/s²     | m      | s       | m/s       |
| Astrodynamics (km/s) | km³/s² | km  | s  | km/s      |
| Canonical   | 1 (DU³/TU²) | DU  | TU      | DU/TU     |

Mixing units (e.g., mu in km³/s² with r in meters) silently produces incorrect velocities —
the crate performs no unit checking or conversion.

Common values:
- Earth GM: 398600.4418 km³/s² = 3.986004418e14 m³/s²
- Solar GM: 1.327124400e11 km³/s²
- Earth canonical DU = 6378.137 km, TU ≈ 806.81 s

## Important Constraints

- `mu` must be positive and finite
- `r1` and `r2` must have non-zero magnitude
- `tof` must be positive
- 180° transfers (vectors exactly anti-parallel) are singular → `LambertError::SingularTransfer`
- Multi-revolution solutions only exist above a minimum TOF; the minimum grows with `nrev`
- `nrev=0` always has a solution (for valid geometry and positive tof)

## Common Mistakes

1. **Mixing units** — use a consistent set throughout (don't mix km and m)
2. **Wrong direction for retrograde** — retrograde uses the long arc; if tof is sized for the
   short arc it will be too small and the geometry changes
3. **Assuming uniqueness for multi-rev** — two solutions exist per revolution count above T_min;
   `lambert` returns the long-period one; there is no short-period entry point yet
4. **180° transfers** — detect with `r1 · r2 / (|r1| · |r2|) ≈ −1` and handle separately

## Finding the Minimum Time of Flight for Multi-Revolution Transfers

No closed-form formula exists for T_min(nrev). To find the minimum time for a given `nrev`,
binary-search on `tof` until `LambertError::NoSolution` transitions to `Ok`:

```rust
// Rough approach: halve tof until NoSolution, then bisect upward
let mut lo = 0.0_f64;
let mut hi = tof_estimate;
for _ in 0..60 {
    let mid = (lo + hi) / 2.0;
    match lambert(mu, r1, r2, mid, nrev, dir) {
        Ok(_) => hi = mid,
        Err(LambertError::NoSolution) => lo = mid,
        _ => break,
    }
}
// `hi` is now approximately T_min; `lambert(..., hi, nrev, dir)` is the minimum-energy solution
```

## Dependencies

This crate has no runtime dependencies — pure Rust `f64` arithmetic only.

Optional feature `gooding-ffi` enables C FFI tests against the original C reference
implementation (requires a C compiler at build time; not needed for library use).

## Algorithm Reference

See [docs/gooding-lambert.md](gooding-lambert.md) for a detailed
description of the four-layer algorithm (tlamb → xlamb → vlamb → lambert), Gooding's
dimensionless orbit parameter x, the Householder iteration, and the velocity recovery scheme.
