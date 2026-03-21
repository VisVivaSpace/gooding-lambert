# gooding-lambert

[![crates.io](https://img.shields.io/crates/v/gooding-lambert.svg)](https://crates.io/crates/gooding-lambert)
[![docs.rs](https://docs.rs/gooding-lambert/badge.svg)](https://docs.rs/gooding-lambert)
[![license](https://img.shields.io/crates/l/gooding-lambert.svg)](LICENSE)

A Rust implementation of Gooding's method for solving Lambert's problem: given two position
vectors and a time of flight, find the connecting orbit's initial and final velocity vectors.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
gooding-lambert = "0.1.1"
```

## Quick Start

```rust
use gooding_lambert::{lambert, Direction};

let mu = 398600.4418_f64;          // Earth GM (km³/s²)
let r1 = [6678.0, 0.0, 0.0];      // LEO departure (km)
let r2 = [0.0, 42164.0, 0.0];     // GEO arrival (km)
let tof = 5.0 * 3600.0;           // 5 hours (seconds)

let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde).unwrap();
println!("Departure v: {:?}", sol.v1);  // km/s
println!("Arrival v:   {:?}", sol.v2);  // km/s
```

## Why Gooding?

Lambert's problem --- finding the orbit connecting two positions in a given time --- is central
to astrodynamics: it underpins interplanetary trajectory design (pork-chop plots), orbit
determination from position observations, rendezvous targeting, and conjunction screening.
Gooding's 1990 method parameterizes the problem through a single dimensionless variable that
spans all conic sections (elliptic, parabolic, hyperbolic) without case-splitting, and solves
the resulting time equation using Householder iteration with analytically computed first, second,
and third derivatives --- no finite differences, no numerical differentiation. Three iterations
reliably yield 13-digit accuracy. The method handles multi-revolution solutions (finding the
minimum-time boundary iteratively, then solving for both long-period and short-period branches)
and degrades gracefully near the parabolic limit via a dedicated power series, making it the
most robust classical Lambert solver available.

## Documentation

- [API Documentation](https://docs.rs/gooding-lambert)
- [Algorithm Details](docs/gooding-lambert.md)
- [LLM Quick-Start Guide](docs/llm-context.md)
- [Gooding (1990) Paper](docs/gooding-paper.pdf)

## Development

This project was co-developed with [Claude](https://claude.ai), an AI assistant by Anthropic.

## License

MIT License — see [LICENSE](LICENSE) for details.
