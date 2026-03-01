# gooding-lambert

A Rust implementation of Gooding's method for solving Lambert's problem: given two position
vectors and a time of flight, find the connecting orbit's initial and final velocity vectors.

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
gooding-lambert = "0.1.0"
```

## Quick Start

```rust
use gooding_lambert::lambert;

// Earth gravitational parameter (m³/s²)
let mu = 3.986004418e14_f64;

// Departure position (m)
let r1 = [7_000_000.0, 0.0, 0.0];
// Arrival position (m)
let r2 = [0.0, 7_000_000.0, 0.0];
// Time of flight (s)
let tof = 3600.0;

let solution = lambert(mu, r1, r2, tof, 0).unwrap();
println!("Departure v: {:?}", solution.v1);
println!("Arrival v:   {:?}", solution.v2);
```

## Documentation

- [API Documentation](https://docs.rs/gooding-lambert)
- [Algorithm Details](docs/algorithms/gooding-lambert.md)

## Development

This project was co-developed with [Claude](https://claude.ai), an AI assistant by Anthropic.

## License

MIT License — see [LICENSE](LICENSE) for details.
