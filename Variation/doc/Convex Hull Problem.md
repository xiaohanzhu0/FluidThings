


# Variational mesh generation with anisotropic metric conformity

---

## 1) Mapping, metric, and generalized energy

We seek a smooth map $x:\hat\Omega\to\Omega\subset\mathbb{R}^d$ from a convex computational box $\hat\Omega=[0,1]^d$ with coordinates $s=(s_1,\dots,s_d)$ into a (possibly non-convex) physical domain with coordinates $x=(x_1,\dots,x_d)$. A given SPD metric field $M(x)\in\mathbb{R}^{d\times d}$ prescribes anisotropic resolution. With computational spacings $\sigma=(\sigma_1,\dots,\sigma_d)$, we write a functional kernel, or **generalized energy**
$$
\mathcal{J}[x;\sigma]\;=\;\int_{\hat\Omega}\Phi\!\big(M(x(s)),\,\nabla_s x(s);\,\sigma\big)\,ds,
$$
for a convex-in-$\nabla_s x$ density $\Phi$ (prototype choices include misfit-squared, linear-in-metric stretch, and optional orthogonality penalties). Dirichlet and natural boundary conditions encode the target boundary map and tangential smoothness. For example, to enforce unit mesh size in metric-induced space, we can let $\Phi = \sum_\alpha ( M_{ij}\, \frac{\partial x_{i}}{\partial {s_\alpha}} \frac{\partial x_{j}}{\partial {s_\alpha}} - 1)^2$

---

## 2) Euler–Lagrange system as a generalized divergence-form elliptic PDE

The first variation $\delta\mathcal{J}=0$ yields the Euler-Lagrangian equation. For multivariable, multicomponent functions, it is given as:
$$\frac{\partial \mathcal{J}}{\partial x_k} - \sum_\alpha^d \frac{\partial}{\partial s_\alpha} \frac{\partial \mathcal{J}}{\partial x_{k,\alpha}} = 0$$
Which can be expressed as a quasi-linear, vector-valued elliptic system in divergence form:
$$
-\partial_{s_\alpha}\!\Big(A^{ij}_{\alpha\beta}(x,\nabla_s x;\sigma)\,x_{j,\beta}\Big)
\;=\;F^i(x,\nabla_s x;\sigma)
\quad\text{in }\hat\Omega,\qquad i=1,\dots,d,
$$
with boundary conditions of Dirichlet/Neumann type inherited from the variational statement. Here $x_{j,\beta}=\partial x_j/\partial s_\beta$, $A^{ij}_{\alpha\beta}$ is strongly elliptic and symmetric in $(i,\alpha)\leftrightarrow(j,\beta)$ under standard kernel choices, and $F^i$ collects lower-order terms. Different practical kernels simply instantiate different $(A,F)$ while preserving this structure.

The original paper solves in the fully-expanded form. We currently find it can perform better to solve in the divergence form.

For example, the Euler-Lagrangian of the alternative cost functional is:
$$2\sum_\alpha \sigma_\alpha^2 M_{kj} \frac{\partial^2 x_j}{\partial s_\alpha^2} + \sum_\alpha \left(2 \frac{\partial M_{kj}}{\partial x_p}\frac{\partial x_p}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha} -  \frac{\partial M_{ij}}{\partial x_k} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha}\right) =0$$
Or:
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left(M_{kj} \frac{\partial x_j}{\partial s_\alpha} \right) = \sum_\alpha \frac{\partial M_{ij}}{\partial x_k} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha}$$

---

## 3) Discretization and linearization

A finite-difference discretization produces a nonlinear system
$$
\mathbf{R}(x)=\mathbf{0},
$$
which we solve by Newton/gradient-like steps (we used Jacobian-Free Newton-Krylov) on $x$, or by a Picard/lagged-coefficient linearization:
$$
\mathcal{L}(x^{(n)})\,\delta x^{(n)}=-\mathbf{R}(x^{(n)}),\qquad x^{(n+1)}=x^{(n)}+\omega^{(n)}\delta x^{(n)}.
$$


---

## 4) The convex-hull phenomenon

### 4.1 What I observe computationally
When mapping $x:\hat\Omega\to\Omega$ with $\hat\Omega$ convex and $\Omega$ non-convex (concave bays, internal obstacles like an airfoil), interior grid lines can “shortcut” across concavities or intrude into obstacles. The numerical image $x(\hat\Omega)$ drifts toward $\operatorname{co}\big(x(\partial\hat\Omega)\big)$ — a **convex-hull effect** — even though the geometric target is non-convex.

### 4.2 Continuous mechanisms (assuming we are working on a uniform metric)
- **Harmonic/elliptic smoothing and convexity:** For linear systems such as $-\Delta_s x=0$, the composition with any affine functional $\ell\circ x$ is harmonic; the boundary maximum principle for all such $\ell$ implies $x(\hat\Omega)\subset\operatorname{co}\big(x(\partial\hat\Omega)\big)$. This is the classical convex-hull property for harmonic maps into $\mathbb{R}^d$.
- **General quasi-linear divergence-form systems:** Our EL system is more general (anisotropy, variable coefficients $M(x)$, nonlinearity in $\nabla_s x$). The convex-hull property need not hold a priori, yet the **smoothing** embodied by strong ellipticity still biases solutions toward “central” values determined by boundary data, especially when obstacle/topology information is absent from the functional.


### 4.3 Theorem
Češík, Antonín. 2023. _Convex Hull Property for Elliptic and Parabolic Systems of PDE_. arXiv. https://arxiv.org/abs/2311.16949
![[Screen Shot 2025-09-11 at 8.42.28 AM.png]]

### 4.4 Example: \[0,1\]^2 to L-shaped Map
Consider the harmonic map from unit square $[0,1]^2$ to a L-shaped domain.
We are solving the Laplace equation with following boundary conditions:
$$u(0,y) = \begin{cases} 1-2y, & y<1/2 \\ 0 & y>=1/2 \end{cases}$$
$$u(1,y) = \begin{cases} 1-y, & y<1/2 \\ 1/2 & y>=1/2 \end{cases}$$
$$u(x,0) = 1 ,\quad u(x,1) = x/2$$
The solution is:

$$u(x,y) = 1-y+\frac{xy}{2} + \sum_{n=1}^\infty \frac{c_n \sinh(n\pi(1-x)) + d_n \sinh(n\pi x)}{\sinh(\pi x)} \sin(n\pi y)$$
$$v(x,y) = y+\frac{x}{2}-\frac{xy}{2} + \sum_{n=1}^\infty \frac{c_n \sinh(n\pi(1-x)) + d_n \sinh(n\pi x)}{\sinh(\pi x)} \sin(n\pi y)$$
With
$$c_n = 2 \int_0^1 g_L(y) \sin(n \pi y) dy = - \frac{4 \sin(n\pi/2) }{(n\pi^2)},\quad d_n = \frac{c_n}{2}$$
$$g_L(y) = \begin{cases} -y,\quad y<1/2 \\ y-1, y\geq 1/2 \end{cases}$$
The resulted map is:
![[Convex_hull_exact.png]]
Which represents a mesh that contains overlap of grids. This is strongly prohibited.


### 4.5  Numerical Example: Airfoil Mesh
Boundary curved is interpolated piecewise linearly.
Metric field is obtained as bicubic interpolation from sampled interior datapoints.
![[outputs/20250627_193358/Figure_06.png]]
![[outputs/20250627_193358/Figure_04.png]]

### 4.6 Diagnosis for our kernels
All kernel variants proposed control **metric conformity and smoothness**, but none encode the **feasible image set** $x(\hat\Omega)\subset\Omega\setminus\Omega_{\text{obs}}$ or the **target topology**. Thus, the energy landscape can reward convex-hull shortcuts.

---

## 5) Our solutions

I’m developing four complementary strategies that preserve ellipticity while constraining the image or encoding topology. Each can be combined with the robust “Approximate” kernel.

### A) Inverse formulation $s(x)$ with comparison principles (scalarization on $\Omega$)
We either transform the previous computational-space functional kernel by change of variables.
$$
\hat{\mathcal{J}}[s;\sigma]\;=\;\int_{\Omega}\Phi\!\big(K(s(x)),\,\nabla_x s(x);\,\sigma\big)|J|^{-1}\, dx,
$$
Or starting form a formulation in physical space:
$$
\hat{\mathcal{J}}[s;\sigma]\;=\;\int_{\Omega}\Phi\!\big(K(s(x)),\,\nabla_x s(x);\,\sigma\big)\, dx,
$$
**PDE:** For $\alpha=1,\dots,d$ solve
$$
-\nabla_x\!\cdot\!\big(|J|^{-1}K(x)\,\nabla_x s_\alpha(x)\big)=0\quad\text{in }\Omega\setminus\Omega_{\text{obs}},
$$
$$
-\nabla_x\!\cdot\!\big(K(x)\,\nabla_x s_\alpha(x)\big)=0\quad\text{in }\Omega\setminus\Omega_{\text{obs}},
$$
with mixed Dirichlet/Neumann data: for a rectangle-like parameterization take $s_1=0,1$ on left/right arcs and natural data on the others; similarly for $s_2$ on bottom/top. Set **constant values on obstacle boundaries** (e.g., $s_1=c_1$, $s_2=c_2$ per component) to ensure wrapping around holes.

**Choice of $K$:** $K(x)\approx M(x)^{-1}$ aligns level sets with metric anisotropy. Uniform ellipticity follows from $\lambda_{\min}(K)\ge\kappa_0>0$.

**Linear solve:** formulating in physical space solves a linear equation. The posterior analysis shows no obvious loss of mesh quality, which is a great advantage.

**Continuous guarantees:** With $K$ symmetric SPD and bounded, each $s_\alpha$ satisfies a strong maximum principle and comparison principle. Obstacles therefore repel level sets; they cannot “cut across” forbidden regions.

This approach has easy obstacle handling. But it solves a slightly different problem, requires solving $d$ scalar PDEs on $\Omega$, which is represented as a curved structured grid. So far, we use finite element on a Q4 mesh.

![[s(x)grid.png]]


---
### B) Two-stages Mapping
$X:$ Physical domain (arbitrarily shaped)
$T:$ Parametrized domain $[0,1]^2$
$S:$ Computational domain $[0,1]^2$

$X\to T$ (or $T\to X$): purely harmonic, only depends on domain geometry
$T \to S$ (or $S\to T$): purely metric conforming, only depends on metric

### Harmonic mapping:
$$\nabla_x^2 t = 0$$
Which is equivalent with Thompson's equation:
$$a x_{k,11} - 2b x_{k,12} + cx_{k,22} = 0, \quad k=1,2$$
$$a = x_{1,2}^2 + x_{2,2}^2,\quad b = x_{1,1}x_{1,2} + x_{2,1}x_{2,2}, \quad c = x_{1,1}^2 + x_{2,1}^2$$
Solving to get the map $x(t)$ and its Jacobian $J_{ij}$



### Metric-conforming mapping
Original Euler-Lagrangian:
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left(M_{kj} \frac{\partial x_j}{\partial s_\alpha} \right) = \sum_\alpha \frac{\partial M_{ij}}{\partial x_k} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha}$$
Substitute $x$ by $t$:
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left(M_{kj}(t(x)) \frac{\partial x_j}{\partial t_i} \frac{\partial t_i}{\partial s_\alpha} \right) = \sum_\alpha \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( \frac{\partial x_i}{\partial t_l} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(\frac{\partial x_j}{\partial t_m} \frac{\partial t_m}{\partial s_\alpha} \right)$$
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left(M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} \right) = \sum_\alpha \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( J_{il} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(J_{jm} \frac{\partial t_m}{\partial s_\alpha} \right)$$

We proposed this approach very currently, and is in the process of verification.



### C) Artificial metric shaping near concave bays
Design $M_\epsilon(x)=M(x)+\mu\,\kappa(x)\,P_t(x)$, where $\kappa$ is a boundary-curvature indicator and $P_t$ projects onto boundary tangents, decreasing tangential diffusion relative to inward normals near concave regions. 

It's just a thought, very difficult and delicate to formulate.

---

### D) Penalization kernel and stabilization
Use the **Approximate** kernel (or the Alternate stretch) to keep the highest-order operator uniformly elliptic. Add vanishing regularization (e.g., an $\eta\,\Delta_s x$ with $\eta$ tied to the current misfit level) to guard against loss of ellipticity as the misfit approaches zero.

Add a penalization kernel such as $P = -\sum_{\alpha} \sigma^2_\alpha \log(x_{i,\alpha}x_{i,\alpha})$ to the alternative cost functional. Then the Euler Lagrangian becomes:
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left((M_{kj}-(x_{i,\alpha}x_{i,\alpha}))^{-1}) \frac{\partial x_j}{\partial s_\alpha} \right) = \sum_\alpha \frac{\partial (M_{ij})}{\partial x_k} \frac{\partial x_i}{\partial s_\alpha} \frac{\partial x_j}{\partial s_\alpha} + \dots$$


However, these penalization kernels cause troubles to stay elliptic. 

We tried to add an $\eta\,\Delta_s x$ to keep the highest-order operator uniformly elliptic, but this regularization make the penalization becomes soft and the convex-hull problem comes back. This dilemma cannot be solved.