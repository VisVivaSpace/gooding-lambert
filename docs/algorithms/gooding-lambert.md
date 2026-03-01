# Gooding's Lambert Solver: Algorithm Reference

## Problem Statement

Lambert's problem: given position vectors **r**₁ and **r**₂ and a time of flight Δt, find the
velocity vectors **v**₁ and **v**₂ at each endpoint of the connecting Keplerian arc. It is the
core routine for orbital transfer design, targeting, and orbit determination.

Gooding's 1990 method reduces the problem to a single scalar equation T(x) = T_target, where
x is Gooding's dimensionless orbit parameter encoding the orbit shape relative to the chord
and semi-perimeter of the transfer triangle. T(x) is evaluated using arctan (elliptic) or ln
(hyperbolic); a power series handles the parabolic limit. The equation is solved by Householder
iteration using analytically computed derivatives T′, T″, T‴.

## Notation Cross-Reference

Some presentations of Gooding's method use λ where this code uses `q` (e.g., Izzo 2015).
The relationship is: `q` (here) = λ (Izzo) = √(r₁r₂)·cos(θ/2)/s. The parameter `qsqfm1`
equals c/s, which equals 1 − q² by Gooding's orbit geometry — the variable name records this
algebraic identity.

---

## Algorithm Structure

The implementation has four nested layers:

```
lambert(mu, r1[3], r2[3], tof, nrev, dir)
  → geometry: th, z (orbit plane normal)
  → vlamb(mu, r1_mag, r2_mag, th, nrev, tof)
      → dimensionless params: q, qsqfm1, rho, sig, T
      → xlamb(m, q, qsqfm1, T, nrev) → x
      → tlamb(m, q, qsqfm1, x, n=-1) → velocity recovery quantities
      → 2D velocities [v_radial, v_transverse] at r1 and r2
  → rotate 2D velocities into 3D Cartesian
```

---

## Layer 1: tlamb — Dimensionless Time Equation

`tlamb(m, q, qsqfm1, x, n)` computes the dimensionless time of flight T(x) and up to three
analytic derivatives, given:

- `x` — Gooding's dimensionless orbit parameter; x ∈ (−1, 1] elliptic, x > 1 hyperbolic
- `m` — revolution count (0 = single arc)
- `q` — geometry parameter = √(r₁r₂) · cos(θ/2) / s
- `qsqfm1` — geometry parameter = c/s (= 1 − q² by Gooding's geometry)

where `s = (r₁ + r₂ + c)/2` is the semi-perimeter and `c` is the chord length.

The `n` argument controls which outputs are computed:

| n  | Outputs                                    | Use                          |
|----|--------------------------------------------|------------------------------|
| 0  | T                                          | Evaluation only              |
| 1  | T, T′                                      |                              |
| 2  | T, T′, T″                                 | Householder-2 refinement     |
| 3  | T, T′, T″, T‴                             | Householder-3 (find T_min)   |
| -1 | (qzminx, qzplx, zplqx) — velocity recovery | After converging on x        |

### Two Code Paths

**Series path** (small-x elliptic region, m=0, x≥0, |1−x²| ≤ SW=0.4):
Accumulates a power series in u = 1 − x² until convergence. Used near the parabolic limit
(x ≈ 0) to avoid catastrophic cancellation in the trigonometric formula.

**Direct path** (all other cases):
- **Elliptic** (x ≤ 1): T via arctan
- **Hyperbolic** (x > 1): T via natural log; uses series log expansion when argument is small
  to avoid cancellation

**Derivatives** are computed via analytic recurrences — no finite differences anywhere.

### n=−1: Velocity Recovery Mode

When called with n=−1 after x has converged, `tlamb` returns three auxiliary quantities
`(qzminx, qzplx, zplqx)` built from z = √(c/s + q²x²) and x. These are combined with
`rho` and `sig` in `vlamb` to recover the 2D velocity components directly, without going
through semi-major axis or Lagrange coefficients.

---

## Layer 2: xlamb — Root Finder

`xlamb(m, q, qsqfm1, T_target, nrev)` finds x such that T(x) = T_target.

### Single Revolution (m=0)

1. **Initial guess**: bilinear/biquadratic approximation based on T_target relative to T(x=0)
2. **Refinement**: 3 iterations of Householder order-2:
   ```
   x ← x + (T_target − T) · T′ / (T′² + (T_target − T) · T″/2)
   ```
3. **Convergence check**: |T_target − T_final| < TOL · max(T_target, 1)

### Multi-Revolution (m>0)

1. **Find T_min**: locate xₘ where T′(xₘ) = 0 (minimum-time point on the m-revolution branch)
   using Householder order-3 on T′:
   ```
   xₘ ← xₘ − T′·T″ / (T″² − T′·T‴/2)
   ```
2. **Check feasibility**: if T_target < T_min(m), return `LambertError::NoSolution`
3. **Select solution branch**: two solutions exist above T_min. The public `lambert` function
   always passes a non-negative `nrev`, so it always returns the long-period solution (x close
   to xₘ, slower orbit). The short-period branch (x further from xₘ) is accessible only by
   calling `xlamb` directly with a negative `nrev`.
4. **Refine** with 3 Householder-2 iterations (same formula as single-rev)

The C reference implementation uses `goto` to share the refinement loop between branches.
The Rust port uses `if/else` to select the initial guess, then a single shared refinement
loop — no gotos required.

---

## Layer 3: vlamb — 2D Velocity Solver

`vlamb(mu, r1, r2, th, nrev, tof)` operates in the 2D radial-transverse plane.

### Dimensionless Parameters

```
dr      = r1 − r2
c       = √((r1−r2)² + 4·r1·r2·sin²(θ/2))   (chord length)
s       = (r1 + r2 + c) / 2                    (semi-perimeter)
gms     = √(μ·s/2)                              (velocity scale)
q       = √(r1·r2) · cos(θ/2) / s
qsqfm1  = c / s
rho     = (r1 − r2) / c                         (radial asymmetry; 0 when r1 = r2)
sig     = 4·r1·r2·sin²(θ/2) / c²               (transverse scale factor)
T       = 4·gms·tof / s²                        (dimensionless time target)
```

### Velocity Recovery

After `xlamb` returns x, call `tlamb(m, q, qsqfm1, x, -1)` to get `(qzminx, qzplx, zplqx)`.
These encode Gooding's velocity coefficients at the solution point. The 2D velocities are:

```
v1_radial     = gms · (qzminx − qzplx · rho) / r1
v1_transverse = gms · zplqx · √sig / r1

v2_radial     = −gms · (qzminx + qzplx · rho) / r2
v2_transverse = gms · zplqx · √sig / r2
```

---

## Layer 4: lambert — 3D Wrapper

`lambert(mu, r1, r2, tof, nrev, dir)` is the public entry point.

### Steps

1. **Validate inputs**: mu > 0, tof > 0, finite vectors, non-zero magnitudes.

2. **Transfer angle**: `th = acos(r1·r2 / (|r1|·|r2|))`, clamped to [−1, 1].
   Detect 180° singularity: if cos θ ≤ −1 + 1e-10, return `LambertError::SingularTransfer`.

3. **Orbit plane normal**: z = r1 × r2 (unit vector). If |z| is tiny (near-collinear),
   return `LambertError::SingularTransfer`.

4. **Retrograde handling**: The prograde transfer uses the short arc (θ ∈ (0°, 180°)).
   For retrograde (θ ∈ (180°, 360°)), two things must change together:
   - `th ← 2π − th`   (use the supplementary arc — changes orbit geometry)
   - `z ← −z`         (flip orbit plane normal — changes angular momentum direction)

5. **Solve 2D**: call `vlamb` with the (possibly modified) th and the scalar magnitudes.
   `nrev > 0` selects the long-period multi-revolution solution.

6. **Build frame**: orthonormal basis at each endpoint:
   ```
   x̂ᵢ = rᵢ / |rᵢ|          (radial unit vector)
   ŷᵢ = ẑ × x̂ᵢ             (transverse unit vector)
   ```

7. **Rotate to 3D**: `vᵢ = v_radial · x̂ᵢ + v_transverse · ŷᵢ`

---

## Numerical Properties

- **Convergence**: Householder-2 is quadratically convergent; 3 iterations are sufficient for
  single-rev from the bilinear guess. Householder-3 for T_min is cubically convergent.
- **Series path**: avoids cancellation near x=0 (parabolic limit). Threshold SW=0.4 is
  generous; the series converges rapidly for |u| ≤ 0.4.
- **Tolerance**: TOL=1e-12 relative to T_target. The public interface achieves ~1e-10 velocity
  agreement with independent solvers (cross-validated against ivLam2's algorithm).
- **No finite differences**: all derivatives T′, T″, T‴ are analytic.
- **No heap allocation**: all intermediate state is on the stack.

---

## References

1. **Gooding, R. H.** (1990). "A procedure for the solution of Lambert's orbital
   boundary-value problem." *Celestial Mechanics and Dynamical Astronomy*, **48**(2), 145–165.
   DOI: [10.1007/BF00049511](https://doi.org/10.1007/BF00049511)

2. **Izzo, D.** (2015). "Revisiting Lambert's problem." *Celestial Mechanics and Dynamical
   Astronomy*, **121**(1), 1–15.
   (Introduces the λ/q reparameterization; useful for cross-validation intuition.)

3. **Gooding, R. H.** (1988). "On the solution of Lambert's orbital boundary-value problem."
   Royal Aircraft Establishment Technical Report 88027. (Precursor to the 1990 paper.)

4. **Battin, R. H.** (1999). *An Introduction to the Mathematics and Methods of
   Astrodynamics*, Revised Ed. AIAA. (Background reference for Keplerian orbit geometry.)

5. **Vallado, D. A.** (2013). *Fundamentals of Astrodynamics and Applications*, 4th Ed.
   Microcosm Press. (Chapter 7 for Lambert problem context and historical survey.)
