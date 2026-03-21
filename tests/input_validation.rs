//! Input validation tests for [`lambert()`].
//!
//! Verifies that `lambert()` rejects invalid inputs with appropriate
//! `LambertError` variants (returns `Err`, does not panic).

use gooding_lambert::{lambert, Direction, LambertError, MultiRevPeriod};

// -- Helpers ------------------------------------------------------------------

/// Default valid inputs for lambert() -- used as baseline; tests override one
/// parameter at a time.
const MU: f64 = 1.0;
const R1: [f64; 3] = [1.0, 0.0, 0.0];
const R2: [f64; 3] = [0.0, 1.0, 0.0];
const TOF: f64 = std::f64::consts::FRAC_PI_2; // 90-degree transfer, valid for nrev=0

fn expect_invalid_input(result: Result<gooding_lambert::LambertSolution, LambertError>) {
    match result {
        Err(LambertError::InvalidInput(_)) => {} // expected
        Err(other) => panic!("expected InvalidInput, got: {other:?}"),
        Ok(_) => panic!("expected Err(InvalidInput), got Ok"),
    }
}

// -- NaN position vectors -----------------------------------------------------

#[test]
fn nan_in_r1_x() {
    let r1 = [f64::NAN, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn nan_in_r1_y() {
    let r1 = [1.0, f64::NAN, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn nan_in_r1_z() {
    let r1 = [1.0, 0.0, f64::NAN];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn nan_in_r2_x() {
    let r2 = [f64::NAN, 0.0, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn nan_in_r2_y() {
    let r2 = [0.0, f64::NAN, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn nan_in_r2_z() {
    let r2 = [0.0, 1.0, f64::NAN];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Zero-magnitude position vectors ------------------------------------------

#[test]
fn zero_r1() {
    let r1 = [0.0, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn zero_r2() {
    let r2 = [0.0, 0.0, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn near_zero_r1() {
    let r1 = [1e-15, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn near_zero_r2() {
    let r2 = [0.0, 1e-15, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- NaN time of flight -------------------------------------------------------

#[test]
fn nan_tof() {
    let result = lambert(MU, R1, R2, f64::NAN, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Negative time of flight --------------------------------------------------

#[test]
fn negative_tof() {
    let result = lambert(MU, R1, R2, -1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Zero time of flight ------------------------------------------------------

#[test]
fn zero_tof() {
    let result = lambert(MU, R1, R2, 0.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Infinity values ----------------------------------------------------------

#[test]
fn inf_tof() {
    let result = lambert(
        MU, R1, R2, f64::INFINITY, 0, Direction::Prograde, MultiRevPeriod::LongPeriod,
    );
    expect_invalid_input(result);
}

#[test]
fn neg_inf_tof() {
    let result = lambert(
        MU, R1, R2, f64::NEG_INFINITY, 0, Direction::Prograde, MultiRevPeriod::LongPeriod,
    );
    expect_invalid_input(result);
}

#[test]
fn inf_in_r1() {
    let r1 = [f64::INFINITY, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn neg_inf_in_r1() {
    let r1 = [f64::NEG_INFINITY, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn inf_in_r2() {
    let r2 = [0.0, f64::INFINITY, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn neg_inf_in_r2() {
    let r2 = [0.0, f64::NEG_INFINITY, 0.0];
    let result = lambert(MU, R1, r2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Invalid mu ---------------------------------------------------------------

#[test]
fn nan_mu() {
    let result = lambert(f64::NAN, R1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn zero_mu() {
    let result = lambert(0.0, R1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn negative_mu() {
    let result = lambert(-1.0, R1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

#[test]
fn inf_mu() {
    let result = lambert(
        f64::INFINITY, R1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod,
    );
    expect_invalid_input(result);
}

// -- Retrograde direction with invalid inputs ---------------------------------
// Verify validation runs before direction encoding (sign flip on tof).

#[test]
fn nan_tof_retrograde() {
    let result = lambert(
        MU, R1, R2, f64::NAN, 0, Direction::Retrograde, MultiRevPeriod::LongPeriod,
    );
    expect_invalid_input(result);
}

#[test]
fn zero_r1_retrograde() {
    let r1 = [0.0, 0.0, 0.0];
    let result = lambert(MU, r1, R2, TOF, 0, Direction::Retrograde, MultiRevPeriod::LongPeriod);
    expect_invalid_input(result);
}

// -- Verify valid inputs succeed (sanity check) --------------------------------
// Ensures the test helper constants actually produce a valid solution.

#[test]
fn baseline_valid_inputs() {
    let result = lambert(MU, R1, R2, TOF, 0, Direction::Prograde, MultiRevPeriod::LongPeriod);
    assert!(result.is_ok(), "baseline valid inputs should succeed: {result:?}");
}
