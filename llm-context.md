# LLM Context for gooding-lambert

## What This Crate Does

Solves Lambert's problem: given two position vectors and a time of flight, compute the initial
and final velocity vectors of the connecting Keplerian orbit. Uses Gooding's (1990) method, which
is robust across all orbit types (elliptic, parabolic, hyperbolic), handles multi-revolution
solutions, and avoids the singularities common in other formulations.

## Key Types and Functions

- `lambert(mu, r1, r2, tof, revs)` — primary entry point; returns `Result<LambertSolution, LambertError>`
- `LambertSolution` — contains `v1: [f64; 3]` and `v2: [f64; 3]` (departure and arrival velocity vectors)
- `LambertError` — error variants for no-solution cases (e.g., 180° transfer, time too short for multi-rev)

## Common Usage Patterns

### Single-revolution transfer

```rust
use gooding_lambert::lambert;

let mu = 3.986004418e14_f64; // Earth GM (m³/s²)
let r1 = [7_000_000.0, 0.0, 0.0];
let r2 = [0.0, 7_000_000.0, 0.0];
let tof = 3600.0; // seconds

let sol = lambert(mu, r1, r2, tof, 0).unwrap();
// sol.v1 = departure velocity, sol.v2 = arrival velocity
```

### Multi-revolution transfer

```rust
// revs = 1 means one full loop before arrival
let sol = lambert(mu, r1, r2, tof, 1)?;
```

## Important Constraints

- All inputs in consistent SI units (meters, seconds)
- `mu` must be positive
- `r1` and `r2` must have non-zero magnitude
- 180° transfers (vectors exactly anti-parallel) are singular and return `LambertError::SingularTransfer`
- Multi-revolution solutions only exist above a minimum time of flight

## Common Mistakes to Avoid

1. Mixing units — use SI throughout (meters, seconds, not km or days)
2. Requesting more revolutions than the time of flight allows — check for `LambertError`
3. Assuming a unique solution for multi-rev cases — Gooding's method can return two solutions per revolution count

## Dependencies Context

This crate has no external dependencies. All math is implemented directly from Gooding's
1990 paper using pure Rust `f64` arithmetic.
