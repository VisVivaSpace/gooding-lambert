# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `MultiRevPeriod` enum (`LongPeriod`, `ShortPeriod`) — `lambert()` now accepts
  a period parameter for multi-rev solutions, giving validated access to both
  solution branches without needing the C-style `gooding_lambert()` API
- `Display` and `std::error::Error` implementations on `LambertError` — enables
  `?` with `anyhow`/`eyre`/`Box<dyn Error>`
- Private `mag3`, `dot3`, `cross3` vector helpers replacing inline computations
- Randomized round-trip tests (5 functions, 100+ cases) with convergence-checked
  Kepler propagator for independent physical verification
- Multi-rev round-trip tests for both `LongPeriod` and `ShortPeriod` (nrev 1–2)
- Input validation tests (26 cases) verifying `lambert()` rejects NaN, Inf,
  zero-magnitude, and negative TOF inputs
- Unit tests for `tlamb` (11 tests covering all 4 code paths) and `xlamb`
  (7 tests including multi-rev bifurcation)
- Asymmetric solver failure counting in cross-validation tests

### Fixed

- `tlamb` power-series loop no longer panics after 200 iterations — returns
  `LambertError::ConvergenceFailed` instead
- `tlamb` hyperbolic log-series loop now has a 200-iteration limit (was unbounded)
- `tlamb` returns `Result` instead of silent `(0,0,0,0)` for near-singular inputs
- `nrev` cast from `u32` to `i32` now uses `try_from` with `InvalidInput` error
  instead of wrapping silently

### Changed

- `lambert()` signature adds `period: MultiRevPeriod` parameter (breaking change)
- `gooding_lambert()` doc comments now warn about lack of validation and direct
  users to `lambert()` as the primary API

### Documentation

- All tolerance thresholds commented as FORTRAN-origin with purpose annotations
- `vlamb_c` documented as test-only infrastructure matching original FORTRAN

## [0.1.1] - 2026-03-20

### Fixed

- Retrograde case error inherited from reference C code — corrected root cause
  and cleaned up cascading compensating errors

### Added

- Gooding (1990) paper included in `docs/` for reference
- Cross-validation tests against Gooding's original Fortran via C translation

## [0.1.0] - 2026-03-01

### Added

- `lambert` function implementing Gooding's (1990) method for Lambert's orbital
  boundary-value problem
- `LambertSolution` struct with departure (`v1`) and arrival (`v2`) velocity vectors
- `LambertError` enum: `InvalidInput`, `SingularTransfer`, `NoSolution`, `ConvergenceFailed`
- `Direction` enum for `Prograde` and `Retrograde` transfers
- Multi-revolution transfers via `nrev: u32` parameter; returns the long-period solution
- Optional `gooding-ffi` feature for C reference-implementation cross-validation
- Algorithm reference documentation (`docs/gooding-lambert.md`)
- LLM quick-start guide (`docs/llm-context.md`)

[Unreleased]: https://github.com/VisVivaSpace/gooding-lambert/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/VisVivaSpace/gooding-lambert/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/VisVivaSpace/gooding-lambert/releases/tag/v0.1.0
