# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/nstrange/gooding-lambert/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/nstrange/gooding-lambert/releases/tag/v0.1.0
