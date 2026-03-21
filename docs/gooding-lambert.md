# Gooding's Lambert Solver: Algorithm Reference

## Background and Historical Context

Lambert's problem is one of the oldest and most important boundary-value problems in
celestial mechanics: given two position vectors **r**_1 and **r**_2 in an inertial frame and a
time of flight Delta-t between them, find the Keplerian orbit that connects them --- that is,
determine the velocity vectors **v**_1 and **v**_2 at each endpoint. The problem was first posed
by Johann Heinrich Lambert in 1761, who also proved the key theorem: the transfer time between
two points on a Keplerian orbit depends only on the semi-major axis, the sum of the radii
r_1 + r_2, and the chord length c between the endpoints --- not on eccentricity or any other
orbital element independently.

Lambert's problem is fundamental to astrodynamics. It appears whenever two position-time
observations must be connected by an orbit:

- **Interplanetary trajectory design**: pork-chop plots scan over departure and arrival dates,
  solving Lambert's problem at each grid point to map the Delta-v landscape for mission planning.
- **Orbit determination**: given two position fixes from tracking data (radar, optical, GPS),
  a Lambert solve recovers the connecting orbit as an initial estimate for a differential
  corrector.
- **Rendezvous and proximity operations**: terminal targeting computes the burn required to
  reach a specified position at a specified time.
- **Conjunction analysis**: screening debris encounters requires rapid orbit transfer
  calculations across large catalogs.

### The Lancaster-Blanchard Formulation (1969)

Lancaster and Blanchard introduced a unified parameterization of the time-of-flight equation
that works across all conic sections (elliptic, parabolic, hyperbolic) without case-splitting.
Their key idea was to express the transfer time as a function of a single scalar parameter
related to the geometry of the "transfer triangle" formed by **r**_1, **r**_2, and the chord
between them. However, their formulation had numerical weaknesses near the parabolic limit
and did not fully address multi-revolution solutions.

### Gooding's Contribution (1988/1990)

R. H. Gooding, working at the Royal Aerospace Establishment (RAE) in Farnborough, published a
precursor technical report (RAE TR 88027, 1988) and then a definitive journal paper (1990) that
resolved the practical shortcomings of the Lancaster-Blanchard approach. Gooding's contributions
were:

1. **A single dimensionless parameter x** that encodes orbit type and shape relative to the
   chord-and-semi-perimeter geometry. The parameter x spans (-1, +inf): values in (-1, 1)
   are elliptic, x = 1 is parabolic, and x > 1 is hyperbolic. This gives a unified
   treatment of all conic sections with no branching on orbit type in the iteration.

2. **A robust time-of-flight function T(x)** built from arctan (elliptic regime) and ln
   (hyperbolic regime), with a carefully constructed power series in u = 1 - x^2 near the
   parabolic limit (x near 0) to avoid the catastrophic cancellation that plagues direct
   evaluation of the trigonometric/logarithmic formulas when the argument is small.

3. **Analytic derivatives T', T'', T'''** computed via compact recurrence relations, enabling
   high-order Householder iteration (Halley's method for T_min location, Householder-2 for
   root refinement) without any finite differences. Three iterations suffice for 13-digit
   convergence in single-revolution cases.

4. **Complete multi-revolution handling**: for m >= 1 revolutions, the time-of-flight curve
   T(x) has a minimum T_min. Two solutions (long-period and short-period) exist for any
   T_target > T_min. Gooding finds T_min via Householder-3 iteration on T'(x) = 0, then
   branches to find each solution with carefully chosen initial guesses.

5. **Direct velocity recovery** from the converged x without intermediate computation of
   semi-major axis, eccentricity, or Lagrange coefficients. Three auxiliary quantities
   derived from the orbit geometry suffice to express the radial and transverse velocity
   components at each endpoint.

6. **Fortran-77 implementation** in three compact subroutines (`TLAMB`, `XLAMB`, `VLAMB`)
   appended to the 1990 paper, providing a complete reference implementation.

Gooding's method is widely regarded as the most robust classical Lambert solver. It has been
adopted in numerous astrodynamics libraries and flight software systems. Izzo (2015) later
introduced a reformulation with faster convergence for single-revolution cases, but Gooding's
method remains the standard for reliability, especially in multi-revolution scenarios.

---

## Problem Statement

### Mathematical Formulation

Given:
- A gravitational parameter mu (for the central body; mu = GM)
- Two position vectors **r**_1, **r**_2 in an inertial frame (with magnitudes r_1, r_2)
- A time of flight Delta-t > 0
- A revolution count m >= 0 (number of complete revolutions before arrival)
- A transfer direction (prograde or retrograde, determining which arc to use)

Find: velocity vectors **v**_1 at **r**_1 and **v**_2 at **r**_2 such that a body departing
**r**_1 with velocity **v**_1 under two-body gravity arrives at **r**_2 with velocity **v**_2
after time Delta-t, having completed m full revolutions.

### Geometric Setup

The two position vectors and the gravitational center define a triangle (the "transfer
triangle") with:

- **Transfer angle** theta: the angle between **r**_1 and **r**_2, measured in the transfer
  plane. For prograde motion theta is in (0, pi); for retrograde, the supplementary angle
  2*pi - theta is used.
- **Chord** c: the straight-line distance between the endpoints,
  c = sqrt((r_1 - r_2)^2 + 4*r_1*r_2*sin^2(theta/2)).
- **Semi-perimeter** s = (r_1 + r_2 + c) / 2, the half-perimeter of the triangle formed by the
  two radii and the chord.

Lambert's theorem states that the transfer time depends only on s, c, and the semi-major axis a
(or equivalently, the orbit's energy). Gooding's reformulation replaces the semi-major axis with
the dimensionless parameter x, yielding a universal time equation T(x) that is monotone (for
single-revolution) or has a single minimum (for multi-revolution), making root-finding
straightforward.

### The Dimensionless Parameter x

Gooding's parameter x is related to the orbit geometry through the auxiliary variable
z = sqrt(qsqfm1 + q^2 * x^2), where:

- q = sqrt(r_1 * r_2) * cos(theta/2) / s --- a geometry parameter bounded in [-1, 1],
  encoding the shape of the transfer triangle. It equals Izzo's lambda parameter.
- qsqfm1 = 1 - q^2 = c/s --- the complementary geometry parameter. The variable name
  in the code records this algebraic identity.

The parameter x maps continuously across all orbit types:
- x in (-1, 1): elliptic orbits, with x = 0 at the parabolic limit
- x = 1: rectilinear (degenerate parabolic)
- x > 1: hyperbolic orbits

As x increases from -1 toward +inf, the orbit transitions from a highly elliptic long-period
transfer through the parabolic case and into increasingly hyperbolic trajectories with shorter
transfer times. This monotonic relationship (for m = 0) is what makes the root-finding problem
well-posed.

### Dimensionless Time Equation

Gooding reduces Lambert's problem to solving the scalar equation:

    T(x) = T_target

where T_target = 4 * sqrt(mu * s / 2) * Delta-t / s^2 is the dimensionless time corresponding
to the physical time of flight, and T(x) is evaluated differently depending on the regime:

**Elliptic case** (x <= 1): T is expressed in terms of arctan via the substitution
y = sqrt(1 - x^2), leading to an expression involving arctan(f/g) + m*pi, where f and g are
combinations of x, q, y, and z. The m*pi term accounts for multi-revolution solutions.

**Hyperbolic case** (x > 1): T is expressed in terms of the natural logarithm, using
y = sqrt(x^2 - 1). When the argument to ln is near 1 (small f), a series expansion of
ln(1 + w) is used to avoid cancellation.

**Near-parabolic case** (x near 0, m = 0): A power series in u = 1 - x^2 is accumulated term
by term, avoiding the catastrophic cancellation that occurs when the arctan formula is
evaluated near the parabolic limit. The series converges rapidly for |u| <= 0.4 (the threshold
SW used in the code).

The final expression in all cases has the form:

    T = 2 * (arctan_or_ln_expression / y + b) / u

where b is an auxiliary quantity depending on q, x, and z, and u = 1 - x^2.

---

## Notation Cross-Reference

Some presentations of Gooding's method use different symbols for the same quantities:

| This code    | Izzo (2015) | Gooding (1990) | Definition                                    |
|--------------|-------------|----------------|-----------------------------------------------|
| `q`          | lambda      | q              | sqrt(r_1*r_2) * cos(theta/2) / s              |
| `qsqfm1`    | 1 - lambda^2| 1 - q^2        | c / s (algebraic identity)                    |
| `x`          | x           | x              | Dimensionless orbit parameter                 |
| `m`          | M           | N              | Revolution count                              |
| `T`          | T           | T              | Dimensionless time of flight                  |
| `rho`        | ---         | rho            | (r_1 - r_2) / c (radial asymmetry)            |
| `sig`        | ---         | sigma          | 4*r_1*r_2*sin^2(theta/2) / c^2                |

---

## Algorithm Structure

The implementation has four nested layers, following the structure of Gooding's original
Fortran-77 subroutines:

```
lambert(mu, r1[3], r2[3], tof, nrev, dir)
  -> geometry: th, z (orbit plane normal)
  -> vlamb(mu, r1_mag, r2_mag, th, nrev, tof)
      -> dimensionless params: q, qsqfm1, rho, sig, T
      -> xlamb(m, q, qsqfm1, T, nrev) -> x
      -> tlamb(m, q, qsqfm1, x, n=-1) -> velocity recovery quantities
      -> 2D velocities [v_radial, v_transverse] at r1 and r2
  -> rotate 2D velocities into 3D Cartesian
```

---

## Layer 1: tlamb --- Dimensionless Time Equation

`tlamb(m, q, qsqfm1, x, n)` computes the dimensionless time of flight T(x) and up to three
analytic derivatives, given:

- `x` --- Gooding's dimensionless orbit parameter; x in (-1, 1] elliptic, x > 1 hyperbolic
- `m` --- revolution count (0 = single arc)
- `q` --- geometry parameter = sqrt(r_1*r_2) * cos(theta/2) / s
- `qsqfm1` --- geometry parameter = c/s (= 1 - q^2 by Gooding's geometry)

where `s = (r_1 + r_2 + c)/2` is the semi-perimeter and `c` is the chord length.

The `n` argument controls which outputs are computed:

| n  | Outputs                                    | Use                          |
|----|--------------------------------------------|------------------------------|
| 0  | T                                          | Evaluation only              |
| 1  | T, T'                                      |                              |
| 2  | T, T', T''                                 | Householder-2 refinement     |
| 3  | T, T', T'', T'''                           | Householder-3 (find T_min)   |
| -1 | (qzminx, qzplx, zplqx) --- velocity recovery | After converging on x     |

### Two Code Paths

**Series path** (small-x elliptic region, m=0, x>=0, |1-x^2| <= SW=0.4):
Accumulates a power series in u = 1 - x^2 until convergence. Used near the parabolic limit
(x near 0) to avoid catastrophic cancellation in the trigonometric formula.

**Direct path** (all other cases):
- **Elliptic** (x <= 1): T via arctan
- **Hyperbolic** (x > 1): T via natural log; uses series log expansion when argument is small
  to avoid cancellation

**Derivatives** are computed via analytic recurrences. On the direct path, the recurrences are:

```
T'   = (3*x*T - 4*(a + q*x*qsqfm1)/z) / u
T''  = (3*T + 5*x*T' + 4*(q/z)^3 * qsqfm1) / u
T''' = (8*T' + 7*x*T'' - 12*(q/z)^5 * x * qsqfm1) / u
```

where a = z - q*x (or its numerically stable equivalent when q*x > 0), z = sqrt(qsqfm1 + q^2*x^2),
and u = 1 - x^2. These recurrences are derived by direct differentiation of the time equation
with respect to x --- no finite differences anywhere.

On the series path, the derivatives are accumulated term by term within the same power-series
loop, then transformed via chain-rule identities to convert from derivatives with respect to u
to derivatives with respect to x.

### n=-1: Velocity Recovery Mode

When called with n=-1 after x has converged, `tlamb` returns three auxiliary quantities
`(qzminx, qzplx, zplqx)` built from z = sqrt(qsqfm1 + q^2*x^2) and x:

```
qzminx = q*z - x        (or numerically stable equivalent)
qzplx  = q*z + x        (or numerically stable equivalent)
zplqx  = z + q*x        (or numerically stable equivalent)
```

When q*x > 0, the "direct" forms q*z - x and z - q*x suffer cancellation, so the code
computes them via the identity a*b = (a^2 - c^2) * ... to avoid loss of significance.
These quantities are combined with `rho` and `sig` in `vlamb` to recover the 2D velocity
components directly, without going through semi-major axis or Lagrange coefficients.

---

## Layer 2: xlamb --- Root Finder

`xlamb(m, q, qsqfm1, T_target, nrev)` finds x such that T(x) = T_target.

### Single Revolution (m=0)

For m = 0, T(x) is monotonically decreasing from +inf (as x -> -1) to 0 (as x -> +inf),
so there is exactly one root for any T_target > 0.

1. **Initial guess**: bilinear/biquadratic approximation based on T_target relative to T(x=0).
   The guess distinguishes two regimes:
   - T_target <= T(0): x is on the elliptic side near parabolic; initial guess is a simple
     rational function of the time difference.
   - T_target > T(0): x is in the hyperbolic regime; the guess uses a more elaborate
     biquadratic formula with empirical constants (C0, C1, C2) tuned for rapid convergence.

2. **Refinement**: 3 iterations of Householder order-2 (Halley's method):
   ```
   x <- x + (T_target - T) * T' / (T'^2 + (T_target - T) * T''/2)
   ```

3. **Convergence check**: |T_target - T_final| < TOL * max(T_target, 1)

### Multi-Revolution (m>0)

For m >= 1, the time-of-flight curve T(x) has a minimum at some x_m in (-1, 1). Two solutions
exist for T_target > T_min, one on each side of x_m (the "long-period" and "short-period"
solutions, corresponding to more circular and more eccentric connecting orbits respectively).

1. **Find T_min**: locate x_m where T'(x_m) = 0 (minimum-time point on the m-revolution branch)
   using Householder order-3 on T':
   ```
   x_m <- x_m - T'*T'' / (T''^2 - T'*T'''/2)
   ```
   The initial guess for x_m is a function of m and the geometry parameter thr2 = atan2(qsqfm1, 2q)/pi.

2. **Check feasibility**: if T_target < T_min(m), return `LambertError::NoSolution`

3. **Select solution branch**: two solutions exist above T_min. The public `lambert` function
   always passes a non-negative `nrev`, so it always returns the long-period solution (x close
   to x_m, slower orbit). The short-period branch (x further from x_m) is accessible only by
   calling `xlamb` directly with a negative `nrev`.

4. **Refine** with 3 Householder-2 iterations (same formula as single-rev)

The C reference implementation uses `goto` to share the refinement loop between branches.
The Rust port uses `if/else` to select the initial guess, then a single shared refinement
loop --- no gotos required.

---

## Layer 3: vlamb --- 2D Velocity Solver

`vlamb(mu, r1, r2, th, nrev, tof)` operates in the 2D radial-transverse plane.

### Dimensionless Parameters

```
dr      = r1 - r2
c       = sqrt((r1-r2)^2 + 4*r1*r2*sin^2(theta/2))   (chord length)
s       = (r1 + r2 + c) / 2                            (semi-perimeter)
gms     = sqrt(mu*s/2)                                  (velocity scale)
q       = sqrt(r1*r2) * cos(theta/2) / s
qsqfm1  = c / s
rho     = (r1 - r2) / c                                (radial asymmetry; 0 when r1 = r2)
sig     = 4*r1*r2*sin^2(theta/2) / c^2                 (transverse scale factor)
T       = 4*gms*tof / s^2                               (dimensionless time target)
```

Note that rho and sig satisfy rho^2 + sig = 1, reflecting the Pythagorean decomposition of
the chord into radial and transverse components.

### Velocity Recovery

After `xlamb` returns x, call `tlamb(m, q, qsqfm1, x, -1)` to get `(qzminx, qzplx, zplqx)`.
These encode Gooding's velocity coefficients at the solution point. The 2D velocities are:

```
v1_radial     = gms * (qzminx - qzplx * rho) / r1
v1_transverse = gms * zplqx * sqrt(sig) / r1

v2_radial     = -gms * (qzminx + qzplx * rho) / r2
v2_transverse = gms * zplqx * sqrt(sig) / r2
```

The velocity scale gms = sqrt(mu*s/2) makes these expressions dimensionally consistent: gms
has units of length/time (i.e., velocity), and the remaining factors are dimensionless ratios
of the orbit geometry. The sign structure (minus on v2_radial) follows from the fact that the
radial velocity component reverses sign between departure and arrival for the portion
attributable to the asymmetry rho.

---

## Layer 4: lambert --- 3D Wrapper

`lambert(mu, r1, r2, tof, nrev, dir)` is the public entry point.

### Steps

1. **Validate inputs**: mu > 0, tof > 0, finite vectors, non-zero magnitudes.

2. **Transfer angle**: `th = acos(r1.r2 / (|r1|*|r2|))`, clamped to [-1, 1].
   Detect 180-degree singularity: if cos(theta) <= -1 + 1e-10, return `LambertError::SingularTransfer`.

3. **Orbit plane normal**: z = r1 x r2 (unit vector). If |z| is tiny (near-collinear),
   return `LambertError::SingularTransfer`.

4. **Retrograde handling**: The prograde transfer uses the short arc (theta in (0, pi)).
   For retrograde (theta in (pi, 2*pi)), two things must change together:
   - `th <- 2*pi - th`   (use the supplementary arc --- changes orbit geometry)
   - `z <- -z`           (flip orbit plane normal --- changes angular momentum direction)

5. **Solve 2D**: call `vlamb` with the (possibly modified) th and the scalar magnitudes.
   `nrev > 0` selects the long-period multi-revolution solution.

6. **Build frame**: orthonormal basis at each endpoint:
   ```
   x_hat_i = r_i / |r_i|          (radial unit vector)
   y_hat_i = z_hat x x_hat_i      (transverse unit vector)
   ```

7. **Rotate to 3D**: `v_i = v_radial * x_hat_i + v_transverse * y_hat_i`

---

## Numerical Properties

- **Convergence**: Householder-2 is cubically convergent (it uses both T' and T''); 3 iterations
  from the bilinear/biquadratic initial guess are sufficient for single-revolution cases.
  Gooding's paper demonstrates that this achieves 13-digit accuracy consistently.
  Householder-3 for T_min is quartically convergent.
- **Series path**: avoids cancellation near x=0 (parabolic limit). Threshold SW=0.4 is
  generous; the series converges rapidly for |u| <= 0.4.
- **Tolerance**: TOL=1e-12 relative to T_target. The public interface achieves ~1e-10 velocity
  agreement with independent solvers (cross-validated against ivLam2's algorithm).
- **No finite differences**: all derivatives T', T'', T''' are analytic.
- **No heap allocation**: all intermediate state is on the stack.
- **Numerical guards**: when q*x > 0, the quantities a = z - q*x and b = q*z - x are
  computed via algebraically equivalent forms that avoid subtraction of nearly equal values.
  Similarly, the hyperbolic log path uses a series expansion of ln(1 + w) when the argument
  is small.

---

## Edge Cases and Limitations

- **180-degree transfers**: when **r**_1 and **r**_2 are antiparallel, the transfer plane is
  undefined (infinitely many planes contain both vectors). The solver returns
  `LambertError::SingularTransfer`. This is an inherent singularity of the two-body Lambert
  problem, not a limitation of Gooding's method.
- **Near-180-degree transfers**: the orbit plane normal r1 x r2 becomes ill-conditioned as the
  transfer angle approaches 180 degrees, leading to large errors in the 3D velocity direction
  even though the 2D solution (radial and transverse components) remains accurate.
- **Multi-revolution minimum time**: T_min for a given revolution count m is not available in
  closed form. The solver finds it iteratively (Householder-3 on T' = 0), which adds
  computational cost compared to single-revolution solves.
- **Very large revolution counts**: the initial guess for x_m degrades as m grows, potentially
  requiring more Householder iterations. The code allows up to 16 iterations for the T_min
  search.

---

## References

1. **Gooding, R. H.** (1990). "A procedure for the solution of Lambert's orbital
   boundary-value problem." *Celestial Mechanics and Dynamical Astronomy*, **48**(2), 145--165.
   DOI: [10.1007/BF00049511](https://doi.org/10.1007/BF00049511)
   ([PDF](gooding-paper.pdf))

2. **Izzo, D.** (2015). "Revisiting Lambert's problem." *Celestial Mechanics and Dynamical
   Astronomy*, **121**(1), 1--15.
   (Introduces the lambda/q reparameterization; useful for cross-validation intuition.)

3. **Gooding, R. H.** (1988). "On the solution of Lambert's orbital boundary-value problem."
   Royal Aircraft Establishment Technical Report 88027. (Precursor to the 1990 paper;
   available at [DTIC](https://apps.dtic.mil/sti/tr/pdf/ADA200383.pdf).)

4. **Lancaster, E. R. and Blanchard, R. C.** (1969). "A unified form of Lambert's theorem."
   NASA TN D-5368. (Original unified time-of-flight formulation upon which Gooding built.)

5. **Battin, R. H.** (1999). *An Introduction to the Mathematics and Methods of
   Astrodynamics*, Revised Ed. AIAA. (Background reference for Keplerian orbit geometry.)

6. **Vallado, D. A.** (2013). *Fundamentals of Astrodynamics and Applications*, 4th Ed.
   Microcosm Press. (Chapter 7 for Lambert problem context and historical survey.)
