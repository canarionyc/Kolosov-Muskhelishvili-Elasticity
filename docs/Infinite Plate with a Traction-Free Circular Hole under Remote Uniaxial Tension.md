# Infinite Plate with a Traction-Free Circular Hole under Remote Uniaxial Tension  
## Muskhelishvili Complex Variable Formulation  

### 1. Problem Statement  
Consider an infinite elastic plate with a circular hole of radius \( R \), subjected to uniform uniaxial tension \( \sigma^\infty \) in the \( x \)-direction at infinity. The hole boundary is traction-free. We aim to determine the stress and displacement fields using Muskhelishvili’s complex potential approach.

---

### 2. Complex Potentials for the Problem  

We introduce two analytic functions \( \varphi(z) \) and \( \psi(z) \) of the complex variable \( z = x + iy \). For an infinite region with a circular hole under remote loading, the general form of the potentials is:

\[
\varphi(z) = A z + \sum_{n=1}^\infty a_n z^{-n}, \quad
\psi(z) = B z + \sum_{n=1}^\infty b_n z^{-n}
\]

where the singular terms vanish as \( |z| \to \infty \), and \( A, B, a_n, b_n \) are complex constants determined from boundary conditions.

---

### 3. Remote Loading Condition  

At infinity (\( |z| \to \infty \)), the stress state is uniaxial tension:

\[
\sigma_{xx}^\infty = \sigma^\infty, \quad \sigma_{yy}^\infty = 0, \quad \tau_{xy}^\infty = 0
\]

Using the relations between stresses and complex potentials:

\[
\sigma_{xx} + \sigma_{yy} = 2\left[\varphi'(z) + \overline{\varphi'(z)}\right]
\]

\[
\sigma_{yy} - \sigma_{xx} + 2i\tau_{xy} = 2\left[\bar{z}\varphi''(z) + \psi'(z)\right]
\]

As \( |z| \to \infty \), \( \varphi'(z) \to A \), \( \psi'(z) \to B \). Substituting the remote stresses:

\[
\sigma^\infty = 2(A + \bar{A}) \implies A + \bar{A} = \frac{\sigma^\infty}{2}
\]

\[
-\sigma^\infty = 2B \implies B = -\frac{\sigma^\infty}{2}
\]

Choosing \( A \) real (since rotation of coordinates can be absorbed), we set:

\[
A = \frac{\sigma^\infty}{4}, \quad B = -\frac{\sigma^\infty}{2}
\]

Thus:

\[
\varphi(z) = \frac{\sigma^\infty}{4}z + \sum_{n=1}^\infty a_n z^{-n}, \quad
\psi(z) = -\frac{\sigma^\infty}{2}z + \sum_{n=1}^\infty b_n z^{-n}
\]

---

### 4. Traction-Free Boundary Condition on the Hole  

On the hole boundary \( |z| = R \), the traction vector \( T = 0 \). In complex form, the boundary condition is:

\[
\varphi(t) + t\overline{\varphi'(t)} + \overline{\psi(t)} = \text{constant} \quad (|t|=R)
\]

where \( t = R e^{i\theta} \). We may take the constant as zero (rigid body motion ignored). Substituting the series forms and using \( \bar{t} = R^2/t \) on the boundary:

\[
\frac{\sigma^\infty}{4}t + \sum_{n=1}^\infty a_n t^{-n} + t\left[\frac{\sigma^\infty}{4} - \sum_{n=1}^\infty n a_n t^{-n-1}\right] + \overline{-\frac{\sigma^\infty}{2}t + \sum_{n=1}^\infty b_n t^{-n}} = 0
\]

Noting \( \overline{t} = R^2/t \), the conjugate terms become:

\[
\overline{-\frac{\sigma^\infty}{2}t} = -\frac{\sigma^\infty}{2}\frac{R^2}{t}, \quad
\overline{b_n t^{-n}} = \bar{b}_n t^n R^{-2n}
\]

Matching coefficients of like powers of \( t \) (after multiplying through by \( t \) to clear negative powers), we obtain:

- Coefficient of \( t \): \( \frac{\sigma^\infty}{4} + \frac{\sigma^\infty}{4} = \frac{\sigma^\infty}{2} \) (balanced by constant adjustment).
- Coefficient of \( t^{-1} \): \( a_1 + R^2\left(\frac{\sigma^\infty}{4}\right) - \frac{\sigma^\infty}{2}R^2 + \bar{b}_1 R^{-2} = 0 \).
- Coefficient of \( t^{-n} \) (\( n \ge 2 \)): \( a_n - (n-2)a_{n-2}R^2 + \bar{b}_n R^{-2n} = 0 \).

Given symmetry and remote uniaxial loading, only \( n=2 \) term survives. Solving:

\[
a_1 = \frac{\sigma^\infty R^2}{2}, \quad b_1 = 0, \quad a_2 = 0, \quad b_2 = \frac{\sigma^\infty R^4}{2}
\]

All other \( a_n, b_n = 0 \).

---

### 5. Final Complex Potentials  

\[
\boxed{\varphi(z) = \frac{\sigma^\infty}{4}\left(z + \frac{2R^2}{z}\right)}
\]

\[
\boxed{\psi(z) = -\frac{\sigma^\infty}{2}\left(z - \frac{R^4}{z^3} + \frac{R^2}{z}\right)}
\]

These potentials satisfy all boundary conditions.

---

### 6. Stress Field Computation  

Differentiate \( \varphi(z) \):

\[
\varphi'(z) = \frac{\sigma^\infty}{4}\left(1 - \frac{2R^2}{z^2}\right), \quad
\varphi''(z) = \frac{\sigma^\infty}{2}\cdot\frac{R^2}{z^3}
\]

Also:

\[
\psi'(z) = -\frac{\sigma^\infty}{2}\left(1 + \frac{3R^4}{z^4} - \frac{R^2}{z^2}\right)
\]

In polar coordinates \( z = re^{i\theta} \), the stress components are:

\[
\sigma_{rr} + i\tau_{r\theta} = \varphi'(z) + \overline{\varphi'(z)} - e^{2i\theta}\left[\bar{z}\varphi''(z) + \psi'(z)\right]
\]

After algebraic manipulation:

\[
\boxed{\sigma_{rr} = \frac{\sigma^\infty}{2}\left(1 - \frac{R^2}{r^2}\right) + \frac{\sigma^\infty}{2}\left(1 - \frac{4R^2}{r^2} + \frac{3R^4}{r^4}\right)\cos 2\theta}
\]

\[
\boxed{\sigma_{\theta\theta} = \frac{\sigma^\infty}{2}\left(1 + \frac{R^2}{r^2}\right) - \frac{\sigma^\infty}{2}\left(1 + \frac{3R^4}{r^4}\right)\cos 2\theta}
\]

\[
\boxed{\tau_{r\theta} = -\frac{\sigma^\infty}{2}\left(1 + \frac{2R^2}{r^2} - \frac{3R^4}{r^4}\right)\sin 2\theta}
\]

---

### 7. Displacement Field  

For plane stress, displacements are given by:

\[
2\mu(u + iv) = \kappa\varphi(z) - z\overline{\varphi'(z)} - \overline{\psi(z)}
\]

where \( \mu \) is the shear modulus, \( \kappa = \frac{3-\nu}{1+\nu} \) for plane stress, \( \kappa = 3-4\nu \) for plane strain, and \( \nu \) is Poisson’s ratio.

Substituting the potentials and simplifying:

**Radial displacement:**
\[
\boxed{u_r = \frac{\sigma^\infty R}{8\mu}\left[\left(\frac{r}{R} + \frac{2R}{r}\right)(\kappa+1)\cos\theta + \left(\frac{4R}{r} - \frac{2R^3}{r^3}\right)\cos 3\theta\right]}
\]

**Tangential displacement:**
\[
\boxed{u_\theta = -\frac{\sigma^\infty R}{8\mu}\left[\left(\frac{r}{R} - \frac{2R}{r}\right)(\kappa-1)\sin\theta + \left(\frac{4R}{r} - \frac{2R^3}{r^3}\right)\sin 3\theta\right]}
\]

---

### 8. Key Observations  

1. **Stress Concentration:** At \( \theta = \pi/2 \) or \( 3\pi/2 \) (points perpendicular to loading direction), on the hole boundary \( r = R \):

\[
\sigma_{\theta\theta} = 3\sigma^\infty
\]

The stress concentration factor is \( K_t = 3 \), independent of material properties (purely geometric).

2. **At \( \theta = 0 \) or \( \pi \):**  
\[
\sigma_{\theta\theta} = -\sigma^\infty
\]
(compression).

3. **Displacements:** The hole deforms into an elliptical shape under tension.

---

### 9. Remarks  

- The solution is exact within linear elasticity.
- For plane strain, replace \( \nu \) with \( \nu/(1-\nu) \) in stress-displacement relations.
- The complex potential method elegantly handles boundary conditions via analytic continuation.

---

**References:**  
Muskhelishvili, N. I., *Some Basic Problems of the Mathematical Theory of Elasticity*, 1953.  
England, A. H., *Complex Variable Methods in Elasticity*, 1971.