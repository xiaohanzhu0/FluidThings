1. Double pendulum transient chaos
2. Driven vs damped stationary fixed state
3. Hamiltonian chaos energy

#### Types of motion
Chaos: deterministic aperiodic motion with a sensitivity to initial conditions
Fixed point (0 dimensional): locally linear
Periodic motion (1 dimensional): $f_0, 2f_0, 3f_0,\dots$
Quasi-periodic motion: (2 dimensional) $f_0,f_1$ where $\frac{f_0}{f_1}$ not rational

#### Types of solutions
1. conservative (no attractors) vs dissipative (attractors)
2. variation in parameter causes bifurcation: discontinuous change in behavior with a continuous change in parameter
3. number of dimensions & symmetries


Rational numbers
1. Integer ratio
2. Periodic decimal
3. countable infinite in number
4. dense in I
5. Lebesgue measure 0
6. Haosdorff dimension 0

Irrational numbers
1. Aperiodic k-ary representations
2. uncountable $|\mathbb{R}\backslash \mathbb{Q}| = 2^{\chi_0}$
3. dense on I
4. Full measure




Fractal dimension definition
$D_F = \frac{\ln(N)}{\ln(1/\epsilon)}$
$D_F$ characteristics how many balls $N$ of size $\epsilon$ it takes to cover the object


Hartman-Grobman Theorem


Koch Curve (Snowflake)
Length: $(\frac{4}{3})^{n}$
Number of segments: $4^n$
Epsilon: $\frac{1}{3^n}$
$D_F = \frac{\ln(4)}{\ln(3)}$



Lorentz Attractor
Center direction: tangent to trajectories, 1
Unstable direction: 1
Stable direction: normal to the attractor, Cantor set, 0.09




Tent map $x$ vs logistic map $y$:
$$x_N = \frac{2}{N}\sin^{-1}(y_N^{1/2})$$
$$x_{N+1} = f(x_N) = g(y_N)$$
Topological equivalent: their map function and inverse is continuous 





## Henon Attractor



# Analysis of nonlinear experimental or observational data
1. Priori: 
	1. Is it autonomous (somewhat relates to whether parameters are steady)? 
	2. How parameters are influenced by environment setup?
	3. Noise expectation: independent identically distributed noise, 60Hz noise
2. Time series observations


Time delay embedding, Rwell Taken's Theorem
Starting from scalar time series $x_0,x_1,x_2,\dots$
Form $\mathbf{x}_i = [x_i,x_{i+\tau},x_{i+2\tau},\dots,x_{i+(d-1)\tau}]$ that is $\tau$ integer time delay
A good embedding has no crossings

Example $d=3$, $\tau=1$




Fitting functions to data
1. Should have sane limits, zeros ...
2. Symmetries should agree with the data

Pade approximate?




Dissipative circular map
$$\Theta_{N+1} = \Theta_N + \Omega + bI_{N} + K\sin(\Theta_N)$$
$$I_{N+1} = bI_N +K\sin(\Theta_N)$$
Special case $b=1$: Jacobian = 1






### Spatialtemporal dynamics
Stuff that happen in PDE
Similarities with time dynamics
1. Fixed points
2. Stability analysis 
3. All other dynamical phyenomion


Turbulence
Deterministic S-T dynamics with a large continuous populated range of aperiodic motion in both frequency and wave number


Kuramoto-Sivashivsky equation
$$\frac{\partial \phi}{\partial t} + \phi \frac{\partial \phi}{\partial x} = -\frac{\partial^2 \phi}{\partial x^2} - \frac{\partial^4 \phi}{\partial x^4}$$
Complex Ginzberg-Landau equation
$$A_t = RA + (1+ib)\nabla^2 A - (1+ic)|A|^2 A$$
which models a wave field with slowly varying in space 
$$U = A(x,t)e^{ik_0x} + A^*(x,t)e^{-ik_0x} + \text{higher orders}$$
Time translation invariant: $A\to A e^{i\phi}$
Reflection invariant: $A(x,t)=A(-x,t)$
Spatial translation invariant: $x\to x+l$



$P(\partial_z v_{z}') \approx \exp(-|\partial_z v_{z}'|)$  


Number of degree of freedom $N=\left( \frac{L}{l_k}\right)^3 = \text{Re}^{9/4}$
Kolmogorov length $l_k$
