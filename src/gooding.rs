//! Gooding's (1990) Lambert algorithm: tlamb → xlamb → vlamb → lambert.

use std::f64::consts::PI;

use crate::{Direction, LambertError, LambertSolution};

const TWOPI: f64 = 2.0 * PI;
const SW: f64 = 0.4;
const TOL: f64 = 1e-12;
const C0: f64 = 1.7;
const C1: f64 = 0.5;
const C2: f64 = 0.03;
const C3: f64 = 0.15;
const C41: f64 = 1.0;
const C42: f64 = 0.24;

/// Dimensionless time-of-flight T(x) and analytic derivatives up to order `n`.
///
/// Returns `(t, dt, d2t, d3t)`.
/// - `n == -1`: velocity recovery mode — returns `(0, qzminx, qzplx, zplqx)`.
/// - `n == 0`: T only (dt, d2t, d3t are 0).
/// - `n == 1`: T and dT/dx.
/// - `n == 2`: T, dT/dx, d²T/dx².
/// - `n == 3`: T, dT/dx, d²T/dx², d³T/dx³.
fn tlamb(m: i32, q: f64, qsqfm1: f64, x: f64, n: i32) -> (f64, f64, f64, f64) {
    let m_f64 = m as f64;

    let qsq = q * q;
    let xsq = x * x;
    let u = (1.0 - x) * (1.0 + x); // = 1 - x²

    // Guard: the direct path divides by u. The series path (m=0, x≥0, n≠-1)
    // handles u→0 correctly; all other cases bail to avoid division by near-zero.
    if u.abs() < 1e-10 && (n == -1 || m != 0 || x < 0.0) {
        return (0.0, 0.0, 0.0, 0.0);
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
            if i >= n && (t - told).abs() < 1e-15 {
                break;
            }
            if i > 200 {
                panic!("tlamb: series did not converge after 200 iterations (x={x}, q={q}, qsqfm1={qsqfm1})");
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

        return (t, dt, d2t, d3t);
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
        return (0.0, b, bb_v, aa_v);
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
        loop {
            twoi1 += 2.0;
            term *= fg1sq;
            let told = t_acc;
            t_acc += term / twoi1;
            if (t_acc - told).abs() < 1e-15 {
                break;
            }
        }
        t_acc
    };

    let t = 2.0 * (t / y + b) / u;
    let mut dt = 0.0;
    let mut d2t = 0.0;
    let mut d3t = 0.0;

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

    (t, dt, d2t, d3t)
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
        let (t0, _, _, _) = tlamb(m, q, qsqfm1, 0.0, 0);
        let tdiff = tin - t0;

        let mut x = if tdiff <= 0.0 {
            // x is on the elliptic side near parabolic
            t0 * tdiff / (-4.0 * tin)
        } else {
            let mut x = -tdiff / (tdiff + 4.0);
            let w = x + C0 * (2.0 * (1.0 - thr2)).sqrt();
            if w < 0.0 {
                x -= (-w).sqrt().sqrt().sqrt().sqrt()
                    * (x + (tdiff / (tdiff + 1.5 * t0)).sqrt());
            }
            let w = 4.0 / (4.0 + tdiff);
            x * (1.0 + x * (C1 * w - C2 * x * w.sqrt()))
        };

        // Householder order-2 refinement
        for _ in 0..3 {
            let (t, dt, d2t, _) = tlamb(m, q, qsqfm1, x, 2);
            let t_err = tin - t;
            if dt.abs() > 1e-30 {
                x += t_err * dt / (dt * dt + t_err * d2t / 2.0);
            }
        }

        // Convergence check
        let (t_final, _, _, _) = tlamb(m, q, qsqfm1, x, 0);
        if (tin - t_final).abs() < TOL * tin.max(1.0) {
            return Ok(x);
        }
        return Err(LambertError::ConvergenceFailed);
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
        let (_, dt, d2t, d3t) = tlamb(m, q, qsqfm1, xm, 3);
        d2t_xm = d2t;
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

    let (tmin, _, _, _) = tlamb(m, q, qsqfm1, xm, 0);
    let tdiffm = tin - tmin;

    if tdiffm < 0.0 {
        return Err(LambertError::NoSolution);
    }
    if tdiffm.abs() < 1e-14 {
        return Ok(xm);
    }

    // Two solutions exist; select based on nrev sign
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
        x = x * (1.0
            - (1.0 + m_f64 + C41 * (thr2 - 0.5)) / (1.0 + C3 * m_f64)
                * x
                * (C1 * w + C2 * x * w.sqrt()))
            + xm;
        x
    } else {
        // Short-period solution (x further from xm)
        let (t0, _, _, _) = tlamb(m, q, qsqfm1, 0.0, 0);
        let tdiff0 = t0 - tmin;
        let tdiff = tin - t0;
        if tdiff <= 0.0 {
            xm - (tdiffm / (d2t2 - tdiffm * (d2t2 / tdiff0 - 1.0 / (xm * xm)))).sqrt()
        } else {
            let mut x = -tdiff / (tdiff + 4.0);
            let w = x + C0 * (2.0 * (1.0 - thr2)).sqrt();
            if w < 0.0 {
                x -= (-w).sqrt().sqrt().sqrt().sqrt()
                    * (x + (tdiff / (tdiff + 1.5 * t0)).sqrt());
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
        let (t, dt, d2t, _) = tlamb(m, q, qsqfm1, x, 2);
        let t_err = tin - t;
        if dt.abs() > 1e-30 {
            x += t_err * dt / (dt * dt + t_err * d2t / 2.0);
        }
    }

    // Convergence check
    let (t_final, _, _, _) = tlamb(m, q, qsqfm1, x, 0);
    if (tin - t_final).abs() < TOL * tin.max(1.0) {
        Ok(x)
    } else {
        Err(LambertError::ConvergenceFailed)
    }
}

/// Solve Lambert in the 2D radial-transverse frame.
///
/// `th` is the transfer angle in radians (0 to 2π; > π for retrograde long-arc).
/// `nrev` is the revolution count with sign: > 0 = long-period, ≤ 0 = short-period.
/// Returns `([v_radial, v_transverse]` at r1 and r2.
fn vlamb(
    gm: f64,
    r1: f64,
    r2: f64,
    th: f64,
    nrev: i32,
    tdelt: f64,
) -> Result<([f64; 2], [f64; 2]), LambertError> {
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

    let (rho, sig) = if c > 1e-14 {
        (dr / c, r1r2th / csq)
    } else {
        (0.0, 1.0)
    };

    let t = 4.0 * gms * tdelt / (s * s);

    let x = xlamb(m, q, qsqfm1, t, nrev)?;

    // Recover velocity components via tlamb(n=-1)
    // Returns (0, qzminx, qzplx, zplqx)
    let (_, qzminx, qzplx, zplqx) = tlamb(m, q, qsqfm1, x, -1);

    let v1 = [
        gms * (qzminx - qzplx * rho) / r1,
        gms * zplqx * sig.sqrt() / r1,
    ];
    let v2 = [
        -gms * (qzminx + qzplx * rho) / r2,
        gms * zplqx * sig.sqrt() / r2,
    ];

    Ok((v1, v2))
}

/// Solve Lambert's problem using Gooding's (1990) method.
///
/// # Parameters
///
/// - `mu`: gravitational parameter (m³/s² or km³/s², consistent with position/time units)
/// - `r1`, `r2`: position vectors at departure and arrival
/// - `tof`: time of flight (must be positive)
/// - `nrev`: number of complete revolutions before arrival (0 = single arc).
///   For multi-revolution transfers, returns the long-period solution.
/// - `dir`: [`Direction::Prograde`] or [`Direction::Retrograde`]
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
) -> Result<LambertSolution, LambertError> {
    // Input validation
    if !mu.is_finite() || mu <= 0.0 {
        return Err(LambertError::InvalidInput("mu must be finite and positive"));
    }
    if !tof.is_finite() || tof <= 0.0 {
        return Err(LambertError::InvalidInput("tof must be finite and positive"));
    }
    for &v in r1.iter().chain(r2.iter()) {
        if !v.is_finite() {
            return Err(LambertError::InvalidInput("position vector contains non-finite value"));
        }
    }

    let r1_mag = (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]).sqrt();
    let r2_mag = (r2[0] * r2[0] + r2[1] * r2[1] + r2[2] * r2[2]).sqrt();

    if r1_mag < 1e-10 * (r2_mag.max(1.0)) {
        return Err(LambertError::InvalidInput("r1 has zero or near-zero magnitude"));
    }
    if r2_mag < 1e-10 * (r1_mag.max(1.0)) {
        return Err(LambertError::InvalidInput("r2 has zero or near-zero magnitude"));
    }

    // Transfer angle from [0, π]
    let cos_th = (r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]) / (r1_mag * r2_mag);
    let cos_th = cos_th.clamp(-1.0, 1.0);

    // Singular at exactly 180° (transfer plane undefined)
    if cos_th <= -1.0 + 1e-10 {
        return Err(LambertError::SingularTransfer);
    }

    let mut th = cos_th.acos(); // ∈ [0, π]

    // Transfer plane normal: z = r1 × r2 (prograde) or r2 × r1 (retrograde)
    let mut z = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let zm = (z[0] * z[0] + z[1] * z[1] + z[2] * z[2]).sqrt();

    if zm < 1e-10 * r1_mag * r2_mag {
        return Err(LambertError::SingularTransfer);
    }
    z[0] /= zm;
    z[1] /= zm;
    z[2] /= zm;

    // Retrograde: use supplementary arc (th → 2π - th) and flip normal
    if dir == Direction::Retrograde {
        th = TWOPI - th;
        z[0] = -z[0];
        z[1] = -z[1];
        z[2] = -z[2];
    }

    // Solve 2D; nrev > 0 → long-period, 0 → single arc
    let nrev_c = nrev as i32;
    let (va1, va2) = vlamb(mu, r1_mag, r2_mag, th, nrev_c, tof)?;

    // Build orthonormal frame: x (radial), y (transverse), z (orbit normal)
    let x1 = [r1[0] / r1_mag, r1[1] / r1_mag, r1[2] / r1_mag];
    let x2 = [r2[0] / r2_mag, r2[1] / r2_mag, r2[2] / r2_mag];

    // y = z × x (transverse direction at each endpoint)
    let y1 = [
        z[1] * x1[2] - z[2] * x1[1],
        z[2] * x1[0] - z[0] * x1[2],
        z[0] * x1[1] - z[1] * x1[0],
    ];
    let y2 = [
        z[1] * x2[2] - z[2] * x2[1],
        z[2] * x2[0] - z[0] * x2[2],
        z[0] * x2[1] - z[1] * x2[0],
    ];

    // Rotate 2D [radial, transverse] into 3D
    let v1 = [
        va1[0] * x1[0] + va1[1] * y1[0],
        va1[0] * x1[1] + va1[1] * y1[1],
        va1[0] * x1[2] + va1[1] * y1[2],
    ];
    let v2 = [
        va2[0] * x2[0] + va2[1] * y2[0],
        va2[0] * x2[1] + va2[1] * y2[1],
        va2[0] * x2[2] + va2[1] * y2[2],
    ];

    Ok(LambertSolution { v1, v2 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn mag(v: [f64; 3]) -> f64 {
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }

    fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    // Specific orbital energy: v²/2 - mu/r
    fn energy(r: [f64; 3], v: [f64; 3], mu: f64) -> f64 {
        let v_sq = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
        v_sq / 2.0 - mu / mag(r)
    }

    #[test]
    fn tlamb_at_zero() {
        // At x=0, T should be finite and positive for prograde (m=0, q>0)
        let (t, dt, d2t, _) = tlamb(0, 0.5, 0.5, 0.0, 2);
        assert!(t.is_finite() && t > 0.0, "T(0) should be positive, got {t}");
        assert!(dt.is_finite(), "dT/dx at 0 should be finite, got {dt}");
        assert!(d2t.is_finite(), "d²T/dx² at 0 should be finite, got {d2t}");
    }

    #[test]
    fn tlamb_velocity_mode() {
        // n=-1 returns (0, qzminx, qzplx, zplqx) — all finite
        let (t0, qzminx, qzplx, zplqx) = tlamb(0, 0.5, 0.5, 0.3, -1);
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
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde).unwrap();
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
        let sol = lambert(mu, r1, r2, 2.0, 0, Direction::Prograde).unwrap();
        let h1 = cross(r1, sol.v1);
        let h2 = cross(r2, sol.v2);
        for i in 0..3 {
            // |h| ~ 1; δh ~ |r|·δv ~ 1e-12
            assert!((h1[i] - h2[i]).abs() < 1e-12, "h[{i}]: {:.6e} vs {:.6e}", h1[i], h2[i]);
        }
    }

    #[test]
    fn prograde_positive_hz() {
        // Prograde orbit: hz = (r × v).z > 0 for r1 in +x, r2 in +y
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde).unwrap();
        let h = cross(r1, sol.v1);
        assert!(h[2] > 0.0, "hz should be positive for prograde, got {}", h[2]);
    }

    #[test]
    fn retrograde_negative_hz() {
        // Retrograde orbit: hz < 0
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, 3.0 * PI / 2.0, 0, Direction::Retrograde).unwrap();
        let h = cross(r1, sol.v1);
        assert!(h[2] < 0.0, "hz should be negative for retrograde, got {}", h[2]);
    }

    #[test]
    fn leo_to_geo() {
        // Physical sanity: LEO→GEO departure speed should be 7–12 km/s
        let mu = 398600.4418;
        let r1 = [6678.0, 0.0, 0.0];
        let r2 = [0.0, 42164.0, 0.0];
        let sol = lambert(mu, r1, r2, 5.0 * 3600.0, 0, Direction::Prograde).unwrap();
        let speed = mag(sol.v1);
        assert!((7.0..=12.0).contains(&speed), "LEO departure speed {speed:.2} km/s out of range");
    }

    #[test]
    fn singular_transfer_180_deg() {
        // 180° transfer (collinear with opposite signs) should fail
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [-1.0, 0.0, 0.0];
        assert_eq!(
            lambert(mu, r1, r2, 1.0, 0, Direction::Prograde),
            Err(LambertError::SingularTransfer)
        );
    }

    #[test]
    fn invalid_mu() {
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        assert!(matches!(lambert(0.0, r1, r2, 1.0, 0, Direction::Prograde), Err(LambertError::InvalidInput(_))));
        assert!(matches!(lambert(-1.0, r1, r2, 1.0, 0, Direction::Prograde), Err(LambertError::InvalidInput(_))));
        assert!(matches!(lambert(f64::NAN, r1, r2, 1.0, 0, Direction::Prograde), Err(LambertError::InvalidInput(_))));
    }

    #[test]
    fn invalid_tof() {
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        assert!(matches!(lambert(mu, r1, r2, 0.0, 0, Direction::Prograde), Err(LambertError::InvalidInput(_))));
        assert!(matches!(lambert(mu, r1, r2, -1.0, 0, Direction::Prograde), Err(LambertError::InvalidInput(_))));
    }

    #[test]
    fn hyperbolic_positive_energy() {
        // Very short TOF => hyperbolic transfer
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 1.0, 0.0];
        let sol = lambert(mu, r1, r2, 0.1, 0, Direction::Prograde).unwrap();
        let e = energy(r1, sol.v1, mu);
        assert!(e > 0.0, "Short TOF should give hyperbolic (e>0), got e={e}");
    }

    #[test]
    fn three_dimensional_transfer() {
        // Out-of-plane transfer: v1 must have non-zero z component
        let mu = 1.0;
        let r1 = [1.0, 0.0, 0.0];
        let r2 = [0.0, 0.8, 0.6];
        let sol = lambert(mu, r1, r2, PI / 2.0, 0, Direction::Prograde).unwrap();
        assert!(sol.v1[2].abs() > 1e-6, "Expected 3D velocity, got v1z={}", sol.v1[2]);
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
        let sol0 = lambert(mu, r1, r2, 20.0, 0, Direction::Prograde).unwrap();
        // 1-rev long period
        let sol1 = lambert(mu, r1, r2, 20.0, 1, Direction::Prograde).unwrap();
        // Energy: 1-rev should be on a different (more circular) orbit
        let e0 = energy(r1, sol0.v1, mu);
        let e1 = energy(r1, sol1.v1, mu);
        assert!((e0 - e1).abs() > 1e-6, "0-rev and 1-rev should have different energies");
    }
}
