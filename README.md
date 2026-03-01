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
use gooding_lambert::{lambert, Direction};

let mu = 398600.4418_f64;          // Earth GM (km³/s²)
let r1 = [6678.0, 0.0, 0.0];      // LEO departure (km)
let r2 = [0.0, 42164.0, 0.0];     // GEO arrival (km)
let tof = 5.0 * 3600.0;           // 5 hours (seconds)

let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde).unwrap();
println!("Departure v: {:?}", sol.v1);  // km/s
println!("Arrival v:   {:?}", sol.v2);  // km/s
```

## Documentation

- [API Documentation](https://docs.rs/gooding-lambert)
- [Algorithm Details](docs/algorithms/gooding-lambert.md)
- [LLM Quick-Start Guide](docs/llm-context.md)

## Development

This project was co-developed with [Claude](https://claude.ai), an AI assistant by Anthropic.

## License

MIT License — see [LICENSE](LICENSE) for details.
