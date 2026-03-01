# Gooding's Lambert Solver

## Overview

Lambert's problem: given position vectors **r**₁ and **r**₂ and a time of flight Δt, find the
velocity vectors **v**₁ and **v**₂ at each endpoint of the connecting Keplerian arc. It is the
core routine for orbital transfer design, targeting, and orbit determination.

Gooding's 1990 method solves this via a single universal variable `x` ∈ (−1, 1] for elliptic
orbits (extended beyond for hyperbolic), using a Householder third-order iteration that requires
analytic first, second, and third derivatives of the time equation T(x).

## Mathematical Background

### Geometry

Given r₁ = |**r**₁|, r₂ = |**r**₂|, and the chord length c = |**r**₂ − **r**₁|, define:

```
s  = (r₁ + r₂ + c) / 2          (semi-perimeter)
λ² = 1 − c/s                     (related to transfer angle)
λ  = √(1 − c/s) × sign(sin(Δν)) (signed; negative for "long way" transfers)
```

where Δν is the true anomaly change. λ = 1 corresponds to a parabolic transfer; λ ∈ (−1, 1)
for elliptic; |λ| > 1 for hyperbolic.

### Universal Variable

Gooding defines `x` such that:

- x = 0: parabolic orbit
- x ∈ (−1, 0): elliptic, short-way
- x ∈ (0, 1): elliptic, long-way
- x > 1: hyperbolic

The normalized time of flight is:

```
T(x) = (1/√|1−x²|) × [ arccos(λx + √(1−λ²)·√(1−x²)) − λx·√(1−x²) ]   (elliptic)
```

with a Stumpff-function form for the general (including hyperbolic) case.

### Householder Iteration

To solve T(x) = T_target, apply the Householder third-order method (cubically convergent):

```
xₙ₊₁ = xₙ − f/f' · [ 1 + (f·f'')/(2·f'²) ] / [ 1 + (f·f'')/(f'²) + (f²·f''')/(6·f'³) ]
```

where f = T(x) − T_target and derivatives T', T'', T''' are computed analytically.

### Recovering Velocities

Once x is found, the semi-major axis a and Lagrange coefficients f, g, ḟ, ġ are computed from
x and the geometry, then:

```
v₁ = (r₂ − f·r₁) / g
v₂ = (ġ·r₂ − r₁) / g
```

### Multi-Revolution Solutions

For N complete revolutions before arrival, the minimum time of flight T_min(N) exists. At T_min
there are two coincident solutions; above T_min there are two distinct solutions with the same N.
Gooding's method handles this by finding the T_min(N) point and branching accordingly.

## Implementation Notes

- All derivatives T', T'', T''' are computed analytically — no finite differences
- The continued-fraction expansion is used for T near x = 1 (parabolic limit) to avoid
  catastrophic cancellation
- Near x = ±1 the Stumpff functions must be computed via series, not the trigonometric form
- 180° transfers (λ → ±1) are detected before iteration and returned as `LambertError::SingularTransfer`
- Convergence is typically achieved in 2–5 Householder iterations for elliptic cases

## Performance

- Single-revolution solve: ~1–5 μs typical
- No heap allocation; all intermediate state is on the stack

## References

1. **Gooding, R. H.** (1990). "A procedure for the solution of Lambert's orbital
   boundary-value problem." *Celestial Mechanics and Dynamical Astronomy*, **48**(2), 145–165.

2. **Izzo, D.** (2015). "Revisiting Lambert's problem." *Celestial Mechanics and Dynamical
   Astronomy*, **121**(1), 1–15.

3. **Battin, R. H.** (1999). *An Introduction to the Mathematics and Methods of
   Astrodynamics*, Revised Ed. AIAA. (Chapter 6 for Lagrange coefficients and universal variables)

4. **Vallado, D. A.** (2013). *Fundamentals of Astrodynamics and Applications*, 4th Ed.
   Microcosm Press. (Chapter 7 for Lambert problem context)
