//! Gooding's (1990) Lambert algorithm: tlamb → xlamb → vlamb → lambert.

use std::f64::consts::PI;

use crate::{Direction, LambertError, LambertSolution};

/// Multi-revolution period selection for [`lambert()`].
///
/// When `nrev >= 1`, two distinct solution families exist for the same
/// geometry and time of flight: *long-period* (lower energy, more circular)
/// and *short-period* (higher energy, more eccentric). This enum selects
/// which family to return.
///
/// For `nrev == 0` (single-revolution), there is only one solution and
/// this parameter is ignored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MultiRevPeriod {
    /// Long-period (lower energy) multi-revolution solution.
    ///
    /// This is the solution on the more circular orbit that completes
    /// the requested number of revolutions.
    LongPeriod,
    /// Short-period (higher energy) multi-revolution solution.
    ///
    /// This is the solution on the more eccentric orbit that completes
    /// the requested number of revolutions.
    ShortPeriod,
}

const TWOPI: f64 = 2.0 * PI;
const SW: f64 = 0.4;
const TOL: f64 = 1e-12;
const C0: f64 = 1.7;
const C1: f64 = 0.5;
const C2: f64 = 0.03;
const C3: f64 = 0.15;
const C41: f64 = 1.0;
const C42: f64 = 0.24;

/// Euclidean magnitude of a 3-vector.
fn mag3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Dot product of two 3-vectors.
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dimensionless time-of-flight T(x) and analytic derivatives up to order `n`.
///
/// Returns `(t, dt, d2t, d3t)`.
/// - `n == -1`: velocity recovery mode — returns `(0, qzminx, qzplx, zplqx)`.
/// - `n == 0`: T only (dt, d2t, d3t are 0).
/// - `n == 1`: T and dT/dx.
/// - `n == 2`: T, dT/dx, d²T/dx².
/// - `n == 3`: T, dT/dx, d²T/dx², d³T/dx³.
fn tlamb(
    m: i32,
    q: f64,
    qsqfm1: f64,
    x: f64,
    n: i32,
) -> Result<(f64, f64, f64, f64), LambertError> {
    let m_f64 = m as f64;

    let qsq = q * q;
    let xsq = x * x;
    let u = (1.0 - x) * (1.0 + x); // = 1 - x²

    // Guard: the direct path divides by u. The series path (m=0, x≥0, n≠-1)
    // handles u→0 correctly; all other cases bail to avoid division by near-zero.
    // FORTRAN-origin: 1e-10 near-zero guard for u = 1-x² singularity
    if u.abs() < 1e-10 && (n == -1 || m != 0 || x < 0.0) {
        return Err(LambertError::ConvergenceFailed);
    }

    // Series path: only for small-x elliptic region near parabolic limit.
    // Direct path used for: n=-1, m>0, x<0, or |u|>SW.
    if n != -1 && m == 0 && x >= 0.0 && u.abs() <= SW {
        // Power-series expansion (near parabolic limit)
        let mut u0i = 1.0;
        let mut u1i = 1.0;
        let mut u2i = 1.0;
        let mut u3i = 1.0;
        let mut term = 4.0;
        let mut tq = q * qsqfm1;
        let mut tqsum = if q < 0.5 {
            1.0 - q * qsq
        } else {
            (1.0 / (1.0 + q) + q) * qsqfm1
        };
        let mut ttmold = term / 3.0;
        let mut t = ttmold * tqsum;
        let mut dt = 0.0;
        let mut d2t = 0.0;
        let mut d3t = 0.0;

        let mut i: i32 = 0;
        loop {
            i += 1;
            let p = i as f64;
            u0i *= u;
            if n >= 1 && i > 1 {
                u1i *= u;
            }
            if n >= 2 && i > 2 {
                u2i *= u;
            }
            if n >= 3 && i > 3 {
                u3i *= u;
            }
            term *= (p - 0.5) / p;
            tq *= qsq;
            tqsum += tq;
            let told = t;
            let tterm = term / (2.0 * p + 3.0);
            let tqterm_base = tterm * tqsum;
            t -= u0i * ((1.5 * p + 0.25) * tqterm_base / (p * p - 0.25) - ttmold * tq);
            ttmold = tterm;
            let tqterm = tqterm_base * p;
            if n >= 1 {
                dt += tqterm * u1i;
            }
            if n >= 2 {
                d2t += tqterm * u2i * (p - 1.0);
            }
            if n >= 3 {
                d3t += tqterm * u3i * (p - 1.0) * (p - 2.0);
            }
            // Loop runs at least n times, then until convergence
            // FORTRAN-origin: 1e-15 convergence criterion for power-series accumulation
            if i >= n && (t - told).abs() < 1e-15 {
                break;
            }
            if i > 200 {
                return Err(LambertError::ConvergenceFailed);
            }
        }

        if n >= 3 {
            d3t = 8.0 * x * (1.5 * d2t - xsq * d3t);
        }
        if n >= 2 {
            d2t = 2.0 * (2.0 * xsq * d2t - dt);
        }
        if n >= 1 {
            dt *= -2.0 * x;
        }
        t /= xsq;

        return Ok((t, dt, d2t, d3t));
    }

    // Direct computation (trig for elliptic, log/series for hyperbolic)
    let y = u.abs().sqrt();
    let z = (qsqfm1 + qsq * xsq).sqrt();
    let qx = q * x;

    // Compute a, b (and aa, bb for n=-1 velocity mode)
    let (a, b) = if qx <= 0.0 {
        (z - qx, q * z - x)
    } else {
        let aa = z + qx;
        let bb = q * z + x;
        (qsqfm1 / aa, qsqfm1 * (qsq * u - xsq) / bb)
    };

    if n == -1 {
        // Velocity recovery mode: return (unused, qzminx, qzplx, zplqx).
        // aa_v and bb_v are the "alternate" pair (opposite of a,b).
        let (aa_v, bb_v) = if qx < 0.0 {
            (qsqfm1 / a, qsqfm1 * (qsq * u - xsq) / b)
        } else {
            // qx >= 0: aa_v = z+qx, bb_v = q*z+x
            (z + qx, q * z + x)
        };
        // qzminx = b, qzplx = bb_v, zplqx = aa_v
        return Ok((0.0, b, bb_v, aa_v));
    }

    // Compute T
    let f = a * y;
    let t = if x <= 1.0 {
        // Elliptic / parabolic
        let g = if qx * u >= 0.0 {
            x * z + q * u
        } else {
            (xsq - qsq * u) / (x * z - q * u)
        };
        m_f64 * PI + f.atan2(g)
    } else if f > SW {
        // Hyperbolic, large argument
        (f + x * z + q * u).ln()
    } else {
        // Hyperbolic, small argument — log series to avoid cancellation
        let g = x * z + q * u;
        let fg1 = f / (g + 1.0);
        let mut term = 2.0 * fg1;
        let fg1sq = fg1 * fg1;
        let mut t_acc = term;
        let mut twoi1 = 1.0;
        let mut i: i32 = 0;
        loop {
            i += 1;
            twoi1 += 2.0;
            term *= fg1sq;
            let told = t_acc;
            t_acc += term / twoi1;
            // FORTRAN-origin: 1e-15 convergence criterion for log-series accumulation
            if (t_acc - told).abs() < 1e-15 {
                break;
            }
            if i > 200 {
                return Err(LambertError::ConvergenceFailed);
            }
        }
        t_acc
    };

    let t = 2.0 * (t / y + b) / u;
    let mut dt = 0.0;
    let mut d2t = 0.0;
    let mut d3t = 0.0;

    // FORTRAN-origin: 1e-30 near-zero guard to skip derivative computation when z is degenerate
    if n >= 1 && z.abs() > 1e-30 {
        let qz = q / z;
        let qz2 = qz * qz;
        let qz3 = qz * qz2;
        dt = (3.0 * x * t - 4.0 * (a + qx * qsqfm1) / z) / u;
        if n >= 2 {
            d2t = (3.0 * t + 5.0 * x * dt + 4.0 * qz3 * qsqfm1) / u;
        }
        if n >= 3 {
            d3t = (8.0 * dt + 7.0 * x * d2t - 12.0 * qz3 * qz2 * x * qsqfm1) / u;
        }
    }

    Ok((t, dt, d2t, d3t))
}

/// Find Gooding's orbit parameter `x` for dimensionless time `tin`.
///
/// `m` is the number of complete revolutions (0 for single arc).
/// `nrev` selects long-period (> 0) or short-period (≤ 0) solution for m > 0.
/// Returns `x` on success.
fn xlamb(m: i32, q: f64, qsqfm1: f64, tin: f64, nrev: i32) -> Result<f64, LambertError> {
    let thr2 = qsqfm1.atan2(2.0 * q) / PI;
    let m_f64 = m as f64;

    if m == 0 {
        // Single-rev: bilinear initial guess
        let (t0, _, _, _) = tlamb(m, q, qsqfm1, 0.0, 0)?;
        let tdiff = tin - t0;

        let mut x = if tdiff <= 0.0 {
            // x is on the elliptic side near parabolic
            t0 * tdiff / (-4.0 * tin)
        } else {
            let mut x = -tdiff / (tdiff + 4.0);
            let w = x + C0 * (2.0 * (1.0 - thr2)).sqrt();
            if w < 0.0 {
                x -= (-w).sqrt().sqrt().sqrt().sqrt() * (x + (tdiff / (tdiff + 1.5 * t0)).sqrt());
            }
            let w = 4.0 / (4.0 + tdiff);
            x * (1.0 + x * (C1 * w - C2 * x * w.sqrt()))
        };

        // Householder order-2 refinement (3 iterations, matching C)
        for _ in 0..3 {
            let (t, dt, d2t, _) = tlamb(m, q, qsqfm1, x, 2)?;
            let t_err = tin - t;
            // FORTRAN-origin: 1e-30 near-zero guard to avoid division by zero in Householder step
            if dt.abs() > 1e-30 {
                x += t_err * dt / (dt * dt + t_err * d2t / 2.0);
            }
        }

        return Ok(x);
    }

    // Multi-rev: find T_min at x_min via Householder order-3 on d²T/dx² = 0
    let mut xm = 1.0 / (1.5 * (m_f64 + 0.5) * PI);
    if thr2 < 0.5 {
        xm *= (2.0 * thr2).sqrt().sqrt().sqrt();
    } else if thr2 > 0.5 {
        xm *= 2.0 - (2.0 - 2.0 * thr2).sqrt().sqrt().sqrt();
    }

    let mut d2t_xm = 0.0;
    for count in 0..16 {
        let (_, dt, d2t, d3t) = tlamb(m, q, qsqfm1, xm, 3)?;
        d2t_xm = d2t;
        // FORTRAN-origin: 1e-30 near-zero guard for d²T convergence check
        if d2t.abs() < 1e-30 {
            break;
        }
        let xm_old = xm;
        xm -= dt * d2t / (d2t * d2t - dt * d3t / 2.0);
        if (xm_old / xm - 1.0).abs() <= TOL {
            break;
        }
        if count == 15 {
            return Err(LambertError::ConvergenceFailed);
        }
    }

    let (tmin, _, _, _) = tlamb(m, q, qsqfm1, xm, 0)?;
    let tdiffm = tin - tmin;

    if tdiffm < 0.0 {
        return Err(LambertError::NoSolution);
    }
    // FORTRAN-origin: 1e-14 near-zero tolerance for T ≈ T_min early exit
    if tdiffm.abs() < 1e-14 {
        return Ok(xm);
    }

    // Two solutions exist; select based on nrev sign
    // FORTRAN-origin: 1e-30 near-zero guard for d²T minimum fallback
    let d2t2 = if d2t_xm.abs() < 1e-30 {
        3.0 * m_f64 * PI // = 6*m*PI/2
    } else {
        d2t_xm / 2.0
    };

    let mut x = if nrev > 0 {
        // Long-period solution (x closer to xm from above)
        let mut x = (tdiffm / (d2t2 + tdiffm / ((1.0 - xm).powi(2)))).sqrt();
        let w = xm + x;
        let w = w * 4.0 / (4.0 + tdiffm) + (1.0 - w).powi(2);
        x = x
            * (1.0
                - (1.0 + m_f64 + C41 * (thr2 - 0.5)) / (1.0 + C3 * m_f64)
                    * x
                    * (C1 * w + C2 * x * w.sqrt()))
            + xm;
        x
    } else {
        // Short-period solution (x further from xm)
        let (t0, _, _, _) = tlamb(m, q, qsqfm1, 0.0, 0)?;
        let tdiff0 = t0 - tmin;
        let tdiff = tin - t0;
        if tdiff <= 0.0 {
            xm - (tdiffm / (d2t2 - tdiffm * (d2t2 / tdiff0 - 1.0 / (xm * xm)))).sqrt()
        } else {
            let mut x = -tdiff / (tdiff + 4.0);
            let w = x + C0 * (2.0 * (1.0 - thr2)).sqrt();
            if w < 0.0 {
                x -= (-w).sqrt().sqrt().sqrt().sqrt() * (x + (tdiff / (tdiff + 1.5 * t0)).sqrt());
            }
            let w = 4.0 / (4.0 + tdiff);
            x * (1.0
                + (1.0 + m_f64 + C42 * (thr2 - 0.5)) / (1.0 + C3 * m_f64)
                    * x
                    * (C1 * w - C2 * x * w.sqrt()))
        }
    };

    // Householder order-2 refinement
    for _ in 0..3 {
        let (t, dt, d2t, _) = tlamb(m, q, qsqfm1, x, 2)?;
        let t_err = tin - t;
        // FORTRAN-origin: 1e-30 near-zero guard to avoid division by zero in Householder step
        if dt.abs() > 1e-30 {
            x += t_err * dt / (dt * dt + t_err * d2t / 2.0);
        }
    }

    // Convergence check
    let (t_final, _, _, _) = tlamb(m, q, qsqfm1, x, 0)?;
    if (tin - t_final).abs() < TOL * tin.max(1.0) {
        Ok(x)
    } else {
        Err(LambertError::ConvergenceFailed)
    }
}

/// Solve Lambert's problem using Gooding's (1990) method.
///
/// This is the primary validated API. For the low-level C-matching interface
/// (useful for cross-validation with the reference implementation), see
/// [`gooding_lambert()`].
///
/// # Parameters
///
/// - `mu`: gravitational parameter (m³/s² or km³/s², consistent with position/time units)
/// - `r1`, `r2`: position vectors at departure and arrival
/// - `tof`: time of flight (must be positive)
/// - `nrev`: number of complete revolutions before arrival (0 = single arc)
/// - `dir`: [`Direction::Prograde`] or [`Direction::Retrograde`]
/// - `period`: [`MultiRevPeriod::LongPeriod`] or [`MultiRevPeriod::ShortPeriod`] —
///   selects which solution family to return for multi-revolution transfers
///   (`nrev >= 1`). Ignored when `nrev == 0`.
///
/// # Errors
///
/// - [`LambertError::InvalidInput`] — non-positive tof, zero position magnitude, or non-finite input
/// - [`LambertError::SingularTransfer`] — transfer angle is exactly 180° (plane undefined)
/// - [`LambertError::NoSolution`] — TOF below minimum for the requested revolution count
/// - [`LambertError::ConvergenceFailed`] — Householder iteration did not converge
pub fn lambert(
    mu: f64,
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    nrev: u32,
    dir: Direction,
    period: MultiRevPeriod,
) -> Result<LambertSolution, LambertError> {
    // --- Input validation (keep in wrapper, not in gooding_lambert) ---
    if !mu.is_finite() || mu <= 0.0 {
        return Err(LambertError::InvalidInput("mu must be finite and positive"));
    }
    if !tof.is_finite() || tof <= 0.0 {
        return Err(LambertError::InvalidInput(
            "tof must be finite and positive",
        ));
    }
    for &v in r1.iter().chain(r2.iter()) {
        if !v.is_finite() {
            return Err(LambertError::InvalidInput(
                "position vector contains non-finite value",
            ));
        }
    }
    let r1_mag = mag3(&r1);
    let r2_mag = mag3(&r2);
    // FORTRAN-origin: 1e-10 relative tolerance for near-zero position magnitude
    if r1_mag < 1e-10 * r2_mag.max(1.0) {
        return Err(LambertError::InvalidInput(
            "r1 has zero or near-zero magnitude",
        ));
    }
    // FORTRAN-origin: 1e-10 relative tolerance for near-zero position magnitude
    if r2_mag < 1e-10 * r1_mag.max(1.0) {
        return Err(LambertError::InvalidInput(
            "r2 has zero or near-zero magnitude",
        ));
    }

    // Check for 180-degree singularity (gooding_lambert uses fallback z=[0,0,1]
    // which gives a "solution" that is physically meaningless — catch it here)
    // FORTRAN-origin: 1e-10 geometric tolerance for near-180° singularity detection
    let cos_th = dot3(&r1, &r2) / (r1_mag * r2_mag);
    if cos_th <= -1.0 + 1e-10 {
        return Err(LambertError::SingularTransfer);
    }
    let cross_z = cross3(&r1, &r2);
    let cross_z_mag = mag3(&cross_z);
    // FORTRAN-origin: 1e-10 geometric tolerance for near-collinear transfer detection
    if cross_z_mag < 1e-10 * r1_mag * r2_mag {
        return Err(LambertError::SingularTransfer);
    }

    // Convert nrev count to signed i32 for gooding_lambert's C convention:
    //   positive nrev -> long-period solution
    //   negative nrev -> short-period solution
    // Safe conversion: u32 revolution counts beyond i32::MAX are physically
    // unreasonable; reject them rather than silently wrapping.
    let nrev_abs: i32 = i32::try_from(nrev).map_err(|_| {
        LambertError::InvalidInput("nrev exceeds maximum supported revolution count")
    })?;
    let nrev_signed = match period {
        MultiRevPeriod::LongPeriod => nrev_abs,
        MultiRevPeriod::ShortPeriod => -nrev_abs,
    };

    // Encode direction as sign of tdelt: positive = prograde, negative = retrograde.
    // gooding_lambert handles retrograde internally (supplementary angle + flipped normal).
    let signed_tdelt = match dir {
        Direction::Prograde => tof,
        Direction::Retrograde => -tof,
    };

    let mut v1 = [0.0_f64; 3];
    let mut v2 = [0.0_f64; 3];
    let code = gooding_lambert(mu, &r1, &r2, nrev_signed, signed_tdelt, &mut v1, &mut v2);
    match code {
        c if c > 0 => Ok(LambertSolution { v1, v2 }),
        0 => Err(LambertError::NoSolution),
        _ => Err(LambertError::ConvergenceFailed),
    }
}

/// Internal helper for [`gooding_lambert`]: solve in the 2D radial-transverse frame
/// using C-matching signed conventions for `nrev` and `tdelt`.
///
/// The argument list intentionally matches the original FORTRAN subroutine VLAMB
/// for cross-validation fidelity with the reference implementation. This function
/// is test-only infrastructure (the public API uses `lambert()` with Rust idioms);
/// do not refactor the argument list.
#[allow(clippy::too_many_arguments)]
fn vlamb_c(
    gm: f64,
    r1: f64,
    r2: f64,
    th: f64,
    nrev: i32,
    tdelt: f64,
    v1: &mut [f64; 2],
    v2: &mut [f64; 2],
) -> i32 {
    let m = nrev.unsigned_abs() as i32;
    let thr2 = th / 2.0;
    let dr = r1 - r2;
    let r1r2th = 4.0 * r1 * r2 * thr2.sin().powi(2);
    let csq = dr * dr + r1r2th;
    let c = csq.sqrt();
    let s = (r1 + r2 + c) / 2.0;
    let gms = (gm * s / 2.0).sqrt();
    let qsqfm1 = c / s;
    let q = (r1 * r2).sqrt() * thr2.cos() / s;

    // FORTRAN-origin: 1e-14 near-zero guard for chord length c (equal-radius case)
    let (rho, sig) = if c > 1e-14 {
        (dr / c, r1r2th / csq)
    } else {
        (0.0, 1.0)
    };

    let t = 4.0 * gms * tdelt / (s * s);

    let x = match xlamb(m, q, qsqfm1, t, nrev) {
        Ok(x) => x,
        Err(LambertError::NoSolution) => return 0,
        Err(_) => return -1,
    };

    let (_, qzminx, qzplx, zplqx) = match tlamb(m, q, qsqfm1, x, -1) {
        Ok(v) => v,
        Err(_) => return -1,
    };

    // C evaluation order: compute intermediate, then divide
    v1[1] = gms * zplqx * sig.sqrt();
    v2[1] = v1[1] / r2;
    v1[1] /= r1;
    v1[0] = gms * (qzminx - qzplx * rho) / r1;
    v2[0] = -gms * (qzminx + qzplx * rho) / r2;

    1
}

/// Low-level Lambert solver matching the C/FORTRAN calling convention.
///
/// # Warning: No input validation
///
/// This function performs **no input validation**. It does not check for
/// non-positive `gm`, zero-length position vectors, non-finite values, or
/// 180-degree singular transfers. Invalid inputs produce undefined numerical
/// results (NaN, infinity, or meaningless velocities) without any error signal.
///
/// **Use [`lambert()`] instead** for a validated, idiomatic Rust API. This
/// function is public for two purposes:
///
/// 1. **Cross-validation** with the original C/FORTRAN implementation — the
///    argument list intentionally mirrors the C `lambert()` function, using
///    signed `nrev` and signed `tdelt` with the same conventions.
/// 2. **Advanced use** where callers have already validated inputs and need
///    direct access to the C-matching interface.
///
/// # Calling convention (C/FORTRAN-matching)
///
/// - `nrev`: signed revolution count. `|nrev|` is the number of complete
///   revolutions. The sign selects the solution family for multi-rev cases:
///   positive = long-period, negative = short-period. For `nrev == 0`, the
///   sign is irrelevant.
/// - `tdelt`: signed time of flight. Positive = prograde transfer,
///   negative = retrograde transfer.
///
/// Returns a status code: 1 on success, 0 if no solution exists, -1 on error.
pub fn gooding_lambert(
    gm: f64,
    r1: &[f64; 3],
    r2: &[f64; 3],
    nrev: i32,
    tdelt: f64,
    v1: &mut [f64; 3],
    v2: &mut [f64; 3],
) -> i32 {
    let rad1 = mag3(r1);
    let rad2 = mag3(r2);

    let dot = dot3(r1, r2);
    let th = (dot / (rad1 * rad2)).clamp(-1.0, 1.0).acos();

    let mut va1 = [0.0_f64; 2];
    let mut va2 = [0.0_f64; 2];

    // Retrograde convention (C-matching): negative tdelt requests retrograde.
    // Use supplementary angle so cos(th/2) goes negative (flips q in vlamb_c),
    // pass |tdelt| since time-of-flight is always positive, and flip the orbit
    // normal for 3D reconstruction.
    let retrograde = tdelt < 0.0;
    let (th, tdelt) = if retrograde {
        (TWOPI - th, tdelt.abs())
    } else {
        (th, tdelt)
    };

    let code = vlamb_c(gm, rad1, rad2, th, nrev, tdelt, &mut va1, &mut va2);

    if code > 0 {
        let x1 = [r1[0] / rad1, r1[1] / rad1, r1[2] / rad1];
        let x2 = [r2[0] / rad2, r2[1] / rad2, r2[2] / rad2];

        let mut z = cross3(&x1, &x2);
        let zm = mag3(&z);

        // FORTRAN-origin: 1e-10 near-zero guard for orbit normal magnitude (fallback to z-axis)
        if zm < 1e-10 {
            z = [0.0, 0.0, 1.0];
        } else {
            z[0] /= zm;
            z[1] /= zm;
            z[2] /= zm;
        }

        // Retrograde: flip orbit normal so angular momentum is antiparallel to r1 × r2.
        if retrograde {
            z[0] = -z[0];
            z[1] = -z[1];
            z[2] = -z[2];
        }

        let y1 = cross3(&z, &x1);
        let y2 = cross3(&z, &x2);

        v1[0] = va1[0] * x1[0] + va1[1] * y1[0];
        v1[1] = va1[0] * x1[1] + va1[1] * y1[1];
        v1[2] = va1[0] * x1[2] + va1[1] * y1[2];
        v2[0] = va2[0] * x2[0] + va2[1] * y2[0];
        v2[1] = va2[0] * x2[1] + va2[1] * y2[1];
        v2[2] = va2[0] * x2[2] + va2[1] * y2[2];
    }

    code
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn mag(v: [f64; 3]) -> f64 {
        mag3(&v)
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        cross3(&a, &b)
    }

    // Specific orbital energy: v²/2 - mu/r
    fn energy(r: [f64; 3], v: [f64; 3], mu: f64) -> f64 {
        let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        v_sq / 2.0 - mu / mag(r)
    }

    #[test]
    fn tlamb_at_zero() {
        // At x=0, T should be finite and positive for prograde (m=0, q>0)
        let (t, dt, d2t, _) = tlamb(0, 0.5, 0.5, 0.0, 2).unwrap();
        assert!(t.is_finite() && t > 0.0, "T(0) should be positive, got {t}");
        assert!(dt.is_finite(), "dT/dx at 0 should be finite, got {dt}");
        assert!(d2t.is_finite(), "d²T/dx² at 0 should be finite, got {d2t}");
    }

    #[test]
    fn tlamb_velocity_mode() {
        // n=-1 returns (0, qzminx, qzplx, zplqx) — all finite
        let (t0, qzminx, qzplx, zplqx) = tlamb(0, 0.5, 0.5, 0.3, -1).unwrap();
        assert_eq!(t0, 0.0);
        assert!(qzminx.is_finite());
        assert!(qzplx.is_finite());
        assert!(zplqx.is_finite());
    }

    #[test]
    fn energy_conserved_prograde() {
        // Energy must be equal at both endpoints (same orbit)
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let e1 = energy(r1, sol.v1, mu);
        let e2 = energy(r2, sol.v2, mu);
        // Lambert converges to |δT| < 1e-12; energy error ~ v·δv ~ 1e-12
        assert!((e1 - e2).abs() < 1e-12, "energy: {e1} vs {e2}");
    }

    #[test]
    fn angular_momentum_conserved() {
        // h = r × v must be the same at both endpoints
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.5, 0.0];
        let sol = lambert(mu, r1, r2, 2.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let h1 = cross(r1, sol.v1);
        let h2 = cross(r2, sol.v2);
        for i in 0..3 {
            // |h| ~ 1; δh ~ |r|·δv ~ 1e-12
            assert!(
                (h1[i] - h2[i]).abs() < 1e-12,
                "h[{i}]: {:.6e} vs {:.6e}",
                h1[i],
                h2[i]
            );
        }
    }

    #[test]
    fn prograde_positive_hz() {
        // Prograde orbit: hz = (r × v).z > 0 for r1 in +x, r2 in +y
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let h = cross(r1, sol.v1);
        assert!(
            h[2] > 0.0,
            "hz should be positive for prograde, got {}",
            h[2]
        );
    }

    #[test]
    fn retrograde_negative_hz() {
        // Retrograde orbit: hz = dot(r1×r2, r1×v1) < 0
        // Same geometric angle as prograde (90°), but opposite direction of travel.
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Retrograde, MultiRevPeriod::LongPeriod).unwrap();
        let h = cross(r1, sol.v1);
        assert!(
            h[2] < 0.0,
            "hz should be negative for retrograde, got {}",
            h[2]
        );
    }

    #[test]
    fn leo_to_geo() {
        // Physical sanity: LEO→GEO departure speed should be 7–12 km/s
        let mu = 398600.4418;
        let r1 = [6678.0, 0.0, 0.0];
        let r2 = [0.0, 42164.0, 0.0];
        let sol = lambert(mu, r1, r2, 5.0 * 3600.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let speed = mag(sol.v1);
        assert!(
            (7.0..=12.0).contains(&speed),
            "LEO departure speed {speed:.2} km/s out of range"
        );
    }

    #[test]
    fn singular_transfer_180_deg() {
        // 180° transfer (collinear with opposite signs) should fail
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-1.0, 0.0, 0.0];
        assert_eq!(
            lambert(mu, r1, r2, 1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::SingularTransfer)
        );
    }

    #[test]
    fn invalid_mu() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        assert!(matches!(
            lambert(0.0, r1, r2, 1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::InvalidInput(_))
        ));
        assert!(matches!(
            lambert(-1.0, r1, r2, 1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::InvalidInput(_))
        ));
        assert!(matches!(
            lambert(f64::NAN, r1, r2, 1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::InvalidInput(_))
        ));
    }

    #[test]
    fn invalid_tof() {
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        assert!(matches!(
            lambert(mu, r1, r2, 0.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::InvalidInput(_))
        ));
        assert!(matches!(
            lambert(mu, r1, r2, -1.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod),
            Err(LambertError::InvalidInput(_))
        ));
    }

    #[test]
    fn hyperbolic_positive_energy() {
        // Very short TOF => hyperbolic transfer
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, 0.1, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let e = energy(r1, sol.v1, mu);
        assert!(e > 0.0, "Short TOF should give hyperbolic (e>0), got e={e}");
    }

    #[test]
    fn three_dimensional_transfer() {
        // Out-of-plane transfer: v1 must have non-zero z component
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 0.8, 0.6];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        assert!(
            sol.v1[2].abs() > 1e-6,
            "Expected 3D velocity, got v1z={}",
            sol.v1[2]
        );
        // Energy still conserved
        let e1 = energy(r1, sol.v1, mu);
        let e2 = energy(r2, sol.v2, mu);
        assert!((e1 - e2).abs() < 1e-12, "3D energy: {e1} vs {e2}");
    }

    #[test]
    fn multi_rev_longer_tof() {
        // Multi-rev solution requires a larger TOF than single-rev
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-1.0, 0.1, 0.0];
        // Single-rev
        let sol0 = lambert(mu, r1, r2, 20.0, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        // 1-rev long period
        let sol1 = lambert(mu, r1, r2, 20.0, 1, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        // Energy: 1-rev should be on a different (more circular) orbit
        let e0 = energy(r1, sol0.v1, mu);
        let e1 = energy(r1, sol1.v1, mu);
        assert!(
            (e0 - e1).abs() > 1e-6,
            "0-rev and 1-rev should have different energies"
        );
    }

    #[test]
    fn lambert_wraps_gooding_lambert() {
        // Verify lambert() and gooding_lambert() produce identical results
        // for the same physical inputs.
        let mu = 398600.4418;
        let r1 = [6678.0, 0.0, 0.0];
        let r2 = [0.0, 42164.0, 0.0];
        let tof = 5.0 * 3600.0;

        // lambert() via wrapper
        let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();

        // gooding_lambert() directly
        let mut v1 = [0.0_f64; 3];
        let mut v2 = [0.0_f64; 3];
        let code = gooding_lambert(mu, &r1, &r2, 0, tof, &mut v1, &mut v2);
        assert!(code > 0);

        // Must be bitwise identical (same code path, no floating-point divergence)
        assert_eq!(
            sol.v1, v1,
            "v1 mismatch between lambert() and gooding_lambert()"
        );
        assert_eq!(
            sol.v2, v2,
            "v2 mismatch between lambert() and gooding_lambert()"
        );
    }

    #[test]
    fn lambert_retrograde_physics() {
        // Retrograde orbit: verify energy conservation, angular momentum
        // conservation, and angular momentum direction (hz < 0).
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let tof = PI / 2.0;

        let sol = lambert(mu, r1, r2, tof, 0, Direction::Retrograde, MultiRevPeriod::LongPeriod).unwrap();

        // Energy conserved between endpoints
        let e1 = energy(r1, sol.v1, mu);
        let e2 = energy(r2, sol.v2, mu);
        assert!((e1 - e2).abs() < 1e-12, "energy: {e1} vs {e2}");

        // Angular momentum conserved
        let h1 = cross(r1, sol.v1);
        let h2 = cross(r2, sol.v2);
        for i in 0..3 {
            assert!(
                (h1[i] - h2[i]).abs() < 1e-12,
                "h[{i}]: {:.6e} vs {:.6e}",
                h1[i],
                h2[i]
            );
        }

        // h_ref = r1 × r2 points in +z; retrograde means dot(h_ref, h_actual) < 0
        let h_ref = cross(r1, r2);
        let dot_h = h_ref[0] * h1[0] + h_ref[1] * h1[1] + h_ref[2] * h1[2];
        assert!(
            dot_h < 0.0,
            "dot(r1×r2, r1×v1) should be negative for retrograde, got {}",
            dot_h
        );
    }

    #[test]
    fn lambert_wraps_gooding_lambert_multirev() {
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-1.0, 0.1, 0.0];
        let tof = 20.0;

        let sol = lambert(mu, r1, r2, tof, 1, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();

        let mut v1 = [0.0_f64; 3];
        let mut v2 = [0.0_f64; 3];
        let code = gooding_lambert(mu, &r1, &r2, 1, tof, &mut v1, &mut v2);
        assert!(code > 0);

        assert_eq!(sol.v1, v1);
        assert_eq!(sol.v2, v2);
    }

    // ====================================================================
    // GL-13: Unit tests for untested tlamb paths
    // ====================================================================

    #[test]
    fn tlamb_sw_boundary_series_path() {
        // Exercise the power-series expansion for x near parabolic limit.
        // Series path triggers when: n != -1, m == 0, x >= 0, |u| <= SW (0.4).
        // x = 0.8 => u = 1 - 0.64 = 0.36, which is <= SW.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q; // 0.75
        let x = 0.8;

        let (t, dt, d2t, d3t) = tlamb(0, q, qsqfm1, x, 3).unwrap();
        assert!(t.is_finite() && t > 0.0, "T should be positive on series path, got {t}");
        assert!(dt.is_finite(), "dT should be finite on series path, got {dt}");
        assert!(d2t.is_finite(), "d2T should be finite on series path, got {d2t}");
        assert!(d3t.is_finite(), "d3T should be finite on series path, got {d3t}");
    }

    #[test]
    fn tlamb_sw_boundary_series_vs_direct_continuity() {
        // At x just below the SW boundary, the series and direct paths should
        // give similar T values, confirming continuity across the branch.
        // SW = 0.4, so u = SW when x = sqrt(1 - 0.4) ≈ 0.7746.
        // x = 0.77 => u ≈ 0.4071 (direct path, just above SW).
        // x = 0.78 => u ≈ 0.3916 (series path, just below SW).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        let (t_direct, _, _, _) = tlamb(0, q, qsqfm1, 0.77, 0).unwrap();
        let (t_series, _, _, _) = tlamb(0, q, qsqfm1, 0.78, 0).unwrap();

        // T(x) is continuous, so nearby x values should give nearby T values.
        // The difference should be small (both values are O(1)).
        let rel_diff = (t_direct - t_series).abs() / t_direct.max(t_series);
        assert!(
            rel_diff < 0.1,
            "Series and direct paths should be continuous at SW boundary: \
             T(0.77)={t_direct}, T(0.78)={t_series}, rel_diff={rel_diff}"
        );
    }

    #[test]
    fn tlamb_sw_boundary_large_q() {
        // Series path with q >= 0.5 exercises the alternate tqsum formula:
        //   tqsum = (1/(1+q) + q) * qsqfm1
        let q = 0.8;
        let qsqfm1 = 1.0 - q * q; // 0.36
        let x = 0.85; // u = 1 - 0.7225 = 0.2775, well within SW

        let (t, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        assert!(t.is_finite() && t > 0.0, "T should be positive for large-q series, got {t}");
    }

    #[test]
    fn tlamb_hyperbolic_large_f() {
        // x > 1 takes the hyperbolic branch. Large x gives large f,
        // exercising the `f > SW` log-branch: ln(f + x*z + q*u).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 2.0; // u = 1-4 = -3, deeply hyperbolic

        let (t, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        assert!(t.is_finite() && t > 0.0, "T should be positive on hyperbolic branch, got {t}");
    }

    #[test]
    fn tlamb_hyperbolic_with_derivatives() {
        // Hyperbolic branch with all derivatives requested (n=3).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 1.5; // hyperbolic

        let (t, dt, d2t, d3t) = tlamb(0, q, qsqfm1, x, 3).unwrap();
        assert!(t.is_finite(), "T should be finite, got {t}");
        assert!(dt.is_finite(), "dT should be finite, got {dt}");
        assert!(d2t.is_finite(), "d2T should be finite, got {d2t}");
        assert!(d3t.is_finite(), "d3T should be finite, got {d3t}");
        // dT/dx should be negative on the hyperbolic branch (T decreases with x for x>1)
        assert!(dt < 0.0, "dT/dx should be negative for hyperbolic x, got {dt}");
    }

    #[test]
    fn tlamb_hyperbolic_small_f_log_series() {
        // Hyperbolic branch with small f exercises the log-series path
        // (f <= SW, the "small argument" hyperbolic sub-branch).
        // x slightly > 1 gives small y = sqrt(|u|) = sqrt(x²-1), hence small f = a*y.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 1.01; // u = 1-1.0201 = -0.0201, |u| = 0.0201, y ≈ 0.1418

        let (t, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        assert!(t.is_finite() && t > 0.0, "T should be positive on hyp small-f path, got {t}");
    }

    #[test]
    fn tlamb_multi_rev() {
        // m > 0 always uses the direct computation path (never the series).
        // Verify T is finite and increases with m (more revolutions = more time).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 0.3;

        let (t0, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        let (t1, _, _, _) = tlamb(1, q, qsqfm1, x, 0).unwrap();
        let (t2, _, _, _) = tlamb(2, q, qsqfm1, x, 0).unwrap();

        assert!(t0.is_finite() && t0 > 0.0, "T(m=0) should be positive, got {t0}");
        assert!(t1.is_finite() && t1 > 0.0, "T(m=1) should be positive, got {t1}");
        assert!(t2.is_finite() && t2 > 0.0, "T(m=2) should be positive, got {t2}");
        // Each additional revolution adds ~pi to T
        assert!(t1 > t0, "T(m=1)={t1} should exceed T(m=0)={t0}");
        assert!(t2 > t1, "T(m=2)={t2} should exceed T(m=1)={t1}");
    }

    #[test]
    fn tlamb_multi_rev_with_derivatives() {
        // Multi-rev path with all derivatives (n=3) on the direct path.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 0.3;

        let (t, dt, d2t, d3t) = tlamb(1, q, qsqfm1, x, 3).unwrap();
        assert!(t.is_finite(), "T should be finite for m=1, got {t}");
        assert!(dt.is_finite(), "dT should be finite for m=1, got {dt}");
        assert!(d2t.is_finite(), "d2T should be finite for m=1, got {d2t}");
        assert!(d3t.is_finite(), "d3T should be finite for m=1, got {d3t}");
    }

    #[test]
    fn tlamb_u_near_zero_guard_velocity_mode() {
        // x ≈ 1 makes u = 1-x² ≈ 0. With n=-1 (velocity mode), the guard
        // should fire and return Err(ConvergenceFailed).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 1.0 - 1e-12; // u ≈ 2e-12, well below 1e-10

        let result = tlamb(0, q, qsqfm1, x, -1);
        assert!(
            matches!(result, Err(LambertError::ConvergenceFailed)),
            "u-near-zero guard should fire for n=-1, got {result:?}"
        );
    }

    #[test]
    fn tlamb_u_near_zero_guard_multi_rev() {
        // x ≈ 1 with m > 0: the guard should fire because m != 0.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = 1.0 - 1e-12;

        let result = tlamb(1, q, qsqfm1, x, 0);
        assert!(
            matches!(result, Err(LambertError::ConvergenceFailed)),
            "u-near-zero guard should fire for m=1, got {result:?}"
        );
    }

    #[test]
    fn tlamb_u_near_zero_guard_negative_x() {
        // x = -(1 - eps) makes u ≈ 0. With x < 0, the guard should fire.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x = -(1.0 - 1e-12);

        let result = tlamb(0, q, qsqfm1, x, 0);
        assert!(
            matches!(result, Err(LambertError::ConvergenceFailed)),
            "u-near-zero guard should fire for negative x near -1, got {result:?}"
        );
    }

    #[test]
    fn tlamb_u_near_zero_safe_series_path() {
        // x ≈ 1 with m=0, x>=0, n>=0 takes the series path, which handles
        // u→0 correctly. This should succeed, NOT trigger the guard.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        // x close to 1 so u = (1-x)(1+x) < 1e-10, but the series path handles it.
        let x = 1.0 - 1e-12; // u ≈ 2e-12

        // m=0, x>=0, n=0: series path condition is met (n!=-1 && m==0 && x>=0)
        // but |u| <= SW is also met, so series path runs.
        // The guard condition is: |u| < 1e-10 && (n==-1 || m!=0 || x<0)
        // Since n=0, m=0, x>0, the guard does NOT fire. Series path runs.
        let result = tlamb(0, q, qsqfm1, x, 0);
        assert!(
            result.is_ok(),
            "Series path should handle u→0 safely for m=0, x>=0, n=0, got {result:?}"
        );
    }

    // ====================================================================
    // GL-14: Unit tests for xlamb root finder
    // ====================================================================

    #[test]
    fn xlamb_single_rev_converges() {
        // Single-rev (m=0) root finder should converge for a typical geometry.
        // Use the same q,qsqfm1 and compute a target T from a known x.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x_ref = 0.3;

        let (t_target, _, _, _) = tlamb(0, q, qsqfm1, x_ref, 0).unwrap();
        let x_found = xlamb(0, q, qsqfm1, t_target, 0).unwrap();

        assert!(
            (x_found - x_ref).abs() < 1e-8,
            "xlamb should recover x={x_ref}, got {x_found}"
        );
    }

    #[test]
    fn xlamb_single_rev_elliptic_side() {
        // tdiff <= 0 means x is on the elliptic side near parabolic (x close to 0).
        // T(0) > tin forces the elliptic initial guess.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // T(0) for these parameters
        let (t0, _, _, _) = tlamb(0, q, qsqfm1, 0.0, 0).unwrap();
        // Pick a target T slightly above T(0) — x will be slightly positive
        // (elliptic side). Wait: tdiff = tin - t0. For tdiff <= 0 we need tin <= t0.
        let tin = t0 * 0.8; // below T(0), so tdiff < 0 -> elliptic initial guess
        let x = xlamb(0, q, qsqfm1, tin, 0).unwrap();

        assert!(x.is_finite(), "Should converge on elliptic side, got {x}");
        // Verify by evaluating T at the found x
        let (t_check, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        assert!(
            (t_check - tin).abs() < 1e-10 * tin,
            "T(x_found) should match target: T={t_check}, target={tin}"
        );
    }

    #[test]
    fn xlamb_single_rev_hyperbolic_side() {
        // tdiff > 0 means x is on the hyperbolic side (x negative or beyond parabolic).
        // T(0) < tin forces the hyperbolic initial guess path.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        let (t0, _, _, _) = tlamb(0, q, qsqfm1, 0.0, 0).unwrap();
        // Pick a target T above T(0) — tdiff > 0 -> hyperbolic side
        let tin = t0 * 1.5;
        let x = xlamb(0, q, qsqfm1, tin, 0).unwrap();

        assert!(x.is_finite(), "Should converge on hyperbolic side, got {x}");
        let (t_check, _, _, _) = tlamb(0, q, qsqfm1, x, 0).unwrap();
        assert!(
            (t_check - tin).abs() < 1e-10 * tin,
            "T(x_found) should match target: T={t_check}, target={tin}"
        );
    }

    #[test]
    fn xlamb_multi_rev_long_period() {
        // Multi-rev (m=1) with nrev > 0 selects the long-period solution.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;
        let x_ref = 0.3;

        // Compute T for m=1 at a reference x, then ask xlamb to find it.
        let (t_target, _, _, _) = tlamb(1, q, qsqfm1, x_ref, 0).unwrap();
        let x_found = xlamb(1, q, qsqfm1, t_target, 1).unwrap(); // nrev=1 -> long period

        assert!(
            (x_found - x_ref).abs() < 1e-6,
            "xlamb long-period should recover x={x_ref}, got {x_found}"
        );
    }

    #[test]
    fn xlamb_multi_rev_short_period() {
        // Multi-rev (m=1) with nrev <= 0 selects the short-period solution.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // For multi-rev, two solutions exist at the same T. Use a T above T_min.
        // Find T_min first by evaluating at xm (the minimum), then pick T > T_min.
        // We use a known geometry that produces a valid multi-rev solution.
        let x_ref = -0.3; // negative x for short-period

        let (t_target, _, _, _) = tlamb(1, q, qsqfm1, x_ref, 0).unwrap();
        let x_found = xlamb(1, q, qsqfm1, t_target, -1).unwrap(); // nrev=-1 -> short period

        assert!(
            (x_found - x_ref).abs() < 1e-6,
            "xlamb short-period should recover x={x_ref}, got {x_found}"
        );
    }

    #[test]
    fn xlamb_multi_rev_bifurcation_two_solutions() {
        // For m > 0, the same T (above T_min) yields two distinct solutions:
        // long-period (x closer to xm from above) and short-period (x further away).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // Pick a T well above T_min for m=1.
        // First, find T_min by finding xm and evaluating T there.
        // xlamb handles this internally; we just need a T large enough.
        let x_ref_long = 0.3;
        let (t_target, _, _, _) = tlamb(1, q, qsqfm1, x_ref_long, 0).unwrap();

        let x_long = xlamb(1, q, qsqfm1, t_target, 1).unwrap();
        let x_short = xlamb(1, q, qsqfm1, t_target, -1).unwrap();

        // Both should be valid roots
        let (t_long, _, _, _) = tlamb(1, q, qsqfm1, x_long, 0).unwrap();
        let (t_short, _, _, _) = tlamb(1, q, qsqfm1, x_short, 0).unwrap();
        assert!(
            (t_long - t_target).abs() < 1e-8 * t_target,
            "Long-period T should match: {t_long} vs {t_target}"
        );
        assert!(
            (t_short - t_target).abs() < 1e-8 * t_target,
            "Short-period T should match: {t_short} vs {t_target}"
        );

        // The two x values should be distinct
        assert!(
            (x_long - x_short).abs() > 1e-6,
            "Long and short period should give different x: long={x_long}, short={x_short}"
        );
    }

    #[test]
    fn xlamb_multi_rev_no_solution_below_tmin() {
        // For m > 0, if tin < T_min, no solution exists -> NoSolution error.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // Use a very small T that is certainly below T_min for m=1.
        let result = xlamb(1, q, qsqfm1, 0.1, 1);
        assert!(
            matches!(result, Err(LambertError::NoSolution)),
            "Below T_min should return NoSolution, got {result:?}"
        );
    }

    #[test]
    fn xlamb_multi_rev_at_tmin_returns_xm() {
        // When tin is exactly T_min, xlamb should return xm (within tolerance).
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // Manually find T_min: run Householder on d²T/dx² = 0 to find xm,
        // then evaluate T(xm). We replicate xlamb's own logic here.
        let mut xm = 1.0 / (1.5 * 1.5 * PI); // m=1
        for _ in 0..20 {
            let (_, dt, d2t, d3t) = tlamb(1, q, qsqfm1, xm, 3).unwrap();
            if d2t.abs() < 1e-30 {
                break;
            }
            xm -= dt * d2t / (d2t * d2t - dt * d3t / 2.0);
        }
        let (tmin, _, _, _) = tlamb(1, q, qsqfm1, xm, 0).unwrap();

        // Request at T_min: should succeed and return xm
        let x_found = xlamb(1, q, qsqfm1, tmin, 1).unwrap();
        assert!(
            (x_found - xm).abs() < 1e-6,
            "At T_min, xlamb should return xm={xm}, got {x_found}"
        );
    }

    #[test]
    fn xlamb_single_rev_roundtrip_via_lambert() {
        // End-to-end: use lambert() to get a solution, then verify xlamb
        // is exercised by checking that the result is physically consistent.
        // This tests xlamb indirectly through the full solver pipeline.
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 2.0, 0.0];
        let tof = 3.0;

        let sol = lambert(mu, r1, r2, tof, 0, Direction::Prograde, MultiRevPeriod::LongPeriod).unwrap();
        let e1 = energy(r1, sol.v1, mu);
        let e2 = energy(r2, sol.v2, mu);
        assert!(
            (e1 - e2).abs() < 1e-11,
            "Energy conservation verifies xlamb convergence: {e1} vs {e2}"
        );
    }

    #[test]
    fn xlamb_multi_rev_short_period_tdiff_positive() {
        // Short-period branch has a sub-case where tdiff > 0 (tin > T(0)).
        // This exercises the negative-x initial guess within the short-period logic.
        let q = 0.5;
        let qsqfm1 = 1.0 - q * q;

        // T(m=1, x=0) — the "zero reference"
        let (t0, _, _, _) = tlamb(1, q, qsqfm1, 0.0, 0).unwrap();
        // Pick a target T well above T(0) for m=1; this puts the short-period
        // root in the tdiff > 0 sub-branch.
        let tin = t0 * 1.2;
        let result = xlamb(1, q, qsqfm1, tin, -1);

        // Should either converge or return NoSolution — never panic.
        assert!(
            result.is_ok() || matches!(result, Err(LambertError::NoSolution | LambertError::ConvergenceFailed)),
            "Short-period tdiff>0 branch should be handled, got {result:?}"
        );
        // If it converged, verify the root
        if let Ok(x) = result {
            let (t_check, _, _, _) = tlamb(1, q, qsqfm1, x, 0).unwrap();
            assert!(
                (t_check - tin).abs() < 1e-8 * tin,
                "Root should satisfy T(x)=tin: T={t_check}, tin={tin}"
            );
        }
    }
}
