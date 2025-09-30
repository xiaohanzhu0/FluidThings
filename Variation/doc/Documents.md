
## Two-stages Mapping
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
$$2\sum_\alpha \sigma_\alpha^2 \frac{\partial}{\partial s_\alpha} \left(M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} \right) = \sum_\alpha \sigma_\alpha^2 \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( J_{il} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(J_{jm} \frac{\partial t_m}{\partial s_\alpha} \right)$$



### Update:
$J$ is not necessarily symmetric, so by multiplying $MJ$ we lost the symmetric property. This means that the PDE no longer follows the variational principle. Here is a trick to recover its symmetry:

Multiply by another Jacobian:
$$2\sum_\alpha \sigma_\alpha^2 J_{kr} \frac{\partial}{\partial s_\alpha}  \left(M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} \right) = \sum_\alpha J_{kr} \sigma_\alpha^2 \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( J_{il} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(J_{jm} \frac{\partial t_m}{\partial s_\alpha} \right)$$
Use product rule the LHS becomes:
$$2\sum_\alpha \sigma_\alpha^2  \frac{\partial}{\partial s_\alpha}  \left( J_{kr}M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} \right) - \sigma_\alpha^2  \frac{\partial J_{kr}}{\partial s_\alpha}  M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} $$
Move the second term to RHS, then the RHS becomes:
$$\sum_\alpha J_{kr} \sigma_\alpha^2 \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( J_{il} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(J_{jm} \frac{\partial t_m}{\partial s_\alpha} \right)+2\sigma_\alpha^2  \frac{\partial J_{kr}}{\partial s_\alpha}  M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} $$
Utilizing $\partial M_{ij} / \partial x_k = (\partial M_{ij} / \partial x_l) J_{ln} J_{nk}$ and $\frac{\partial J_{kr}}{\partial s_\alpha}=\frac{\partial J_{kr}}{\partial t_q}\frac{\partial t_q}{\partial s_\alpha}$ , the RHS becomes:
$$\sum_\alpha J_{kr} \sigma_\alpha^2 \frac{\partial M_{ij}(x(t))}{\partial x_k} \left( J_{il} \frac{\partial t_l}{\partial s_\alpha}\right)  \left(J_{jm} \frac{\partial t_m}{\partial s_\alpha} \right)+2\sigma_\alpha^2  \frac{\partial J_{kr}}{\partial t_q}   \frac{\partial t_q}{\partial s_\alpha} M_{kj}(t(x)) J_{ji} \frac{\partial t_i}{\partial s_\alpha} $$

Use another product rule such that:
$$\frac{\partial}{\partial t}(J^T M J) = \frac{\partial J^T}{\partial t} M J + J^T \frac{\partial M}{\partial t} J + J^T M\frac{\partial J}{\partial t} = \frac{\partial J^T}{\partial t} M J + J^T \frac{\partial M}{\partial x}J J + J^T M\frac{\partial J}{\partial t}$$
The first and third term yield the same result. So in index notation:
$$\sum_\alpha \sigma_\alpha \frac{\partial}{\partial t_r} (J^T MJ)_{pq} \frac{\partial t_p}{\partial s_\alpha} \frac{\partial t_q}{\partial s_\alpha}$$
Then eventually we get:
$$2\sum_\alpha \sigma_\alpha \frac{\partial}{\partial s_\alpha} \left[ (J^T MJ)_{rq} \frac{\partial t_q}{\partial s_\alpha}\right] = \sum_\alpha \sigma_\alpha \frac{\partial}{\partial t_r} (J^T MJ)_{pq} \frac{\partial t_p}{\partial s_\alpha} \frac{\partial t_q}{\partial s_\alpha}$$
Which is an elliptic equation that follows the variational principle. Compared with the original PDE in $x$ and $s$, we can think of this as the change of basis by transforming the metric tensor $M$ based on $J$.

This can be equivalently interpreted as: transforming the original metric tensor $M$ using the rank-2 covariant tensor transformation rule ($M$ components in $t$ is equal to $J^T MJ$ in $x$).


## Playful Toy Examples
Now we easily can make mixtures of toy examples to explore with.

$X(T)$ on a uniform grid $T$:
![[Lshape_uniform.png]]
This is equivalent to the full map $S(X)$ if the metric field is uniform.

----

$T(S)$ on a uniform $S$, given the metric $J^T M J = \begin{bmatrix} (1+15t_1)^{-2} & 0 \\ 0 & 1 \end{bmatrix}$

![[Lshape_tgrid.png]]

The composition result $X(S)$:
![[Lshape_stretched.png]]
Which is more equally distributed and better minimize $||misfit||_\infty$ than the original approach.

We can see this approach has the potential to offset the misfit caused by the geometric shape of the domain, if we can artificially modify the transformed metric field.









Regarding the $p^{-1}$ penalization in the $x(s)$ formulation: we noted that the penalization modifies the divergence form to  
$\frac{\partial}{\partial s} \big( (1 - p^{-2}) \frac{\partial x}{\partial s} \big) = \text{Source}$.  
This can become singular if $(1 - p^{-2}) \leq 0$. To avoid this, I moved the $p^{-2}$ term to the right-hand side. This prevents the singularity issue and instead applies a softer penalty on shrunken or overlapping grids. In practice, the grid converged under this formulation, but overlapping problems still remained.

---
On the $s(x)$ formulation: since the optimized mesh size satisfies $l_i M_{ij} l_j = 1$, in matrix form we have $\nabla x^T M \nabla x = I$,  which implies $(\nabla x^T M \nabla x)^{-1} = \nabla s^T M^{-1} \nabla s = I$. Thus, if perfect optimization is achievable, we can redefine $p$ as $s_{\alpha,i} M_{ij} s_{\alpha,j}$. Meanwhile, with $ds = |J|^{-1} dx$, there remains a scalar weight of $|J|^{-1}$. This leads to the Euler–Lagrange equation

$$ \frac{\partial}{\partial x_i} ( J^{-1} M_{ij} \frac{\partial s_\alpha}{\partial x_j}) = 0$$

with no source term, since $\frac{\partial L}{\partial s} = 0$.

So far, I have tested a linearized version by assuming $|J|^{-1} \approx \text{Area of element} \approx \sqrt{\det(M)}$, making the Euler–Lagrange equation linear and solvable in one shot. The results differ somewhat in regions with high metric variation or high curvature. This approach is more robust, requires no iteration, and avoids concave overlap issues.






Next step might be testing on the nonlinear case using real J^{-1}. Also, one shows we can transform back to x(s) systems by utilizing the fact that $\xi_x = |J|^{-1} y_\eta,  \xi_y = -|J|^{-1} x_\eta, \eta_x = -|J|^{-1} y_\xi, \eta_y = |J|^{-1} x_\xi$. This is different from the original Euler-Lagrangian as it is based from s(x) formulation, so it avoids the concave overlap problem while still solved on a uniform cartesian mesh:

[http://persson.berkeley.edu/math228b/slides/meshgeneration_slides.pdf](http://persson.berkeley.edu/math228b/slides/meshgeneration_slides.pdf)


$$J[s] = \int_\Omega \sum_{\alpha=1}^d \sigma_\alpha^2 \frac{\partial s_\alpha}{\partial x_i} M_{ij}^{-1} \frac{\partial s_\alpha}{\partial x_j} \,dx $$







### Jacobian-Free Newton Treatment
Even we fully linearize the system using Jacobian-free with Newton-Krylov approximation, it still causes overlaps. I am thinking about this is the inevitable flaw of the alternative functional kernel.

### Penalize small or overlapped mesh grid
Since the alternative cost function is minimized when $p_\alpha$ shrinks to 0, it doesn't punish collapsed mesh grids and can cause overlaps. We add punishment functional kernels.
$$L_{1} = \sum_{\alpha} \sigma^2_\alpha p_\alpha^{-1}$$
$$L_{2} = -\sum_{\alpha} \sigma^2_\alpha \log(p_\alpha)$$
$$L_{3} = \sum_{\alpha} \sigma^2_\alpha |J(x,s)|^{-1}$$
$$L_{3} = -\sum_{\alpha} \sigma^2_\alpha \log(|J(x,s)|)$$
They don't perform well. If there are mesh grids that are too small, they overshoot so much.

![[failed_example.gif]]

Then we want to make the punishment functional kernels finite. Then one of the easiest returns back to the exact/approximate cost functional. For now, it better than alternative at least on the uniform metric tensor field:

![[approximate.png]]

![[alternative.png]]

But as the metric field has greater variance, the approximate one failed to converge. I need to apply blurring convolution kernel about 50 times to get a metric field smooth enough.

IDEA:
I am thinking about a better way to smoothen the metric field:
https://www.sciencedirect.com/science/article/abs/pii/S0168874X09000912
This is more physically correct as it conserves the positive definiteness of metric field and more conformed to a real mesh.





### Mesh generation on C-shaped airfoil

#### Metric tensor smoothing and interpolation
A 1D Gaussian kernel is used to smoothen the raw metric tensor data in both direction in sequence. 

Since the metric data is not sampled uniformly, we used MATLAB's Delaunay triangulation interpolation function `scatteredInterpolant` (which can be improved in the future since the sampled data is structured, but we don't want because the efficiency improvement is not significant so far).

#### Initial mesh generation
Initially, we kept the use of transfinite interpolation as the initial mesh. It couldn't converge well for the case of airfoil. Then we tried to use hyperbolic mesh generation as the initial mesh, which is based on the formulation:
$$x_\xi x_\eta + y_\xi y_\eta = 0,\quad x_\xi y_\eta - y_\xi x_\eta = J^{-1}$$
After linearization, this becomes:
$$x^{n+1}_\xi x_\eta + x_\xi x_\eta^{n+1} + y_\xi^{n+1}y_\eta + y_\xi y_\eta^{n+1}=0,\quad x^{n+1}_\xi y_\eta + x_\xi y_\eta^{n+1} + y_\xi^{n+1}x_\eta + y_\xi x_\eta^{n+1}= J^{-1} + (J^{n+1})^{-1} $$
Or:
$$\mathbf{Aw}_\xi + \mathbf{Bw}_\eta = \mathbf{f}$$
Where:
$$\mathbf{w} = \begin{bmatrix}x^{n+1} \\ y^{n+1}\end{bmatrix}, \quad \mathbf{A} = \begin{bmatrix} x_\eta & y_\eta \\ y_\eta &-x_\eta\end{bmatrix} , \quad \mathbf{B} = \begin{bmatrix} x_\xi & y_\xi \\ -y_\xi &x_\xi\end{bmatrix}, \quad \mathbf{f}=\begin{bmatrix}0 \\ 2J^{-1}\end{bmatrix}$$
Using central differencing in $\xi$ direction and marching in $\eta$ direction using first order biased differencing.


I tried to introduce artificial diffusion by adding second order derivative to make the grids far away more regular, but my current implementation doesn't do well. I instead added "artificially-artificial diffusion" by smoothing each layer as: $$\mathbf{x}_{i,j} = (\mathbf{x}_{i+1,j}+2\mathbf{x}_{i,j}+\mathbf{x}_{i-1,j})/4$$

#### Airfoil result
`Nx1=435; Nx2=142; grade_ratio=1.01; relaxation=0.1`


![[Variation/outputs/20250623_164858/Figure_06.png]]
![[Variation/outputs/20250623_164858/Figure_05.png]]
![[Variation/outputs/20250623_164858/Figure_04.png]]
![[Variation/outputs/20250623_164858/Figure_03.png]]
![[Variation/outputs/20250623_164858/Figure_02.png]]
![[Variation/outputs/20250623_164858/Figure_01.png]]



For now, the convergence condition is very delicate. It only converges to acceptable result under some specific setups. Other results include: 
1. Intruding airfoil heading edge due to large curvature (no enough tangent grid points),
2. Overlap at the normal boundary from the trailing edge (no enough normal grid points),
3. Overlap at the airfoil heading edge (too many tangent grid points)
4. Blows up due to bad aspect ratio at far field (too many normal grid points),





The recurrence relation is:
$$x_{n+1} = x_n + \alpha(x_n-x_{n-1}), \quad \text{where $x_1$ and $x_N$ is given}$$
The general solution is:
$$x_n=\frac{x_N (\alpha^n-\alpha )+x_1 \left(\alpha ^N-\alpha ^n\right)}{\alpha ^N-\alpha }$$





### Validation case: rotating metric tensor
By rotating the system, a diagonal tensor can contain off-diagonal terms:
$$
\begin{array}{cc}
 M_{11} = \cos (\theta ) \left(M_{11} \cos (\theta )-M_{12} \sin (\theta )\right)-\sin (\theta ) \left(M_{12} \cos (\theta )-M_{22} \sin (\theta )\right) \\
 M_{12} = \cos (\theta ) \left(M_{11} \sin (\theta )+M_{12} \cos (\theta )\right)-\sin (\theta ) \left(M_{12} \sin (\theta )+M_{22} \cos (\theta )\right) \\
 M_{22} = \sin (\theta ) \left(M_{11} \sin (\theta )+M_{12} \cos (\theta )\right)+\cos (\theta ) \left(M_{12} \sin (\theta )+M_{22} \cos (\theta )\right) \\
\end{array}
$$
We rotate the domain and the metric tensor by a same angle while keep my physical coordinate fixed, the correct result should be the same as the original case after rotating back. Now the metric tensor would naturally has off-diagonal terms from tensor transformation rules.

### Find $M(x)$ and $\frac{\partial M(x)}{\partial x}$ numerically
For each iteration, we need the values of $M(x)$ and $\frac{\partial M(x)}{\partial x}$ based on $x$ from previous iteration:
We directly using bicubic interpolation to obtain metric tensor value $M$ based on structured grid points $x_1,x_2$.
Compute the gradient $\nabla_{x_1} M$ and $\nabla_{x_2} M$ as:
$$\begin{bmatrix}\frac{\partial x_1}{\partial s_1} & \frac{\partial x_2}{\partial s_1} \\ \frac{\partial x_1}{\partial s_2} & \frac{\partial x_2}{\partial s_2}\end{bmatrix}\begin{bmatrix}\nabla_{x_1}M \\  \nabla_{x_2}M \end{bmatrix}= J_{ij}\begin{bmatrix}\nabla_{x_1}M \\  \nabla_{x_2}M \end{bmatrix}_{ij}= \begin{bmatrix} \frac{\partial M}{\partial s_1} \\ \frac{\partial M}{\partial s_2}\end{bmatrix}_{ij}$$
For every grid points $(i,j)$. $J$ and $\frac{\partial M}{\partial s}$ is computed from finite difference with respect to $s$. $J$ is $2\times 2$, which is very cheap to invert. 
![[2.png]]
![[1.png]]



### Next step: testing on airfoil with multi-block partition
Smoothing by doing convolution in physical domain vs in computational domain
![[example.png]]

Smoothing original sampled data (more consistent, better convergence) vs smoothing each interpolation (more accurate)

Connect work with Marvyn: where to partition
	partition at huge metric spikes?







Based on problem 1: $$M_{11}=40000(1+15x_1)^{-2},\quad M_{22}=40000(1+15x_2)^{-2}$$
$$M_{12}=2000(1+15x_1x_2)^{-2}, \quad\text{or}\quad M_{12}=0$$

![[off_diag.png]]
![[off_diag (1).png]]

![[no_off_diag.png]]
![[no_off_diag (1).png]]

![[p1_residual.png]]




Based on problem 3: $$M_{11}=M_{22}=2000$$
$$M_{12}=1990,\quad \text{or} \quad M_{12}=0$$
With off-diagonal:
![[p3_offdiag.png]]
Without off-diagonal:
![[p3_no_offdiag.png]]
![[p3_residual.png]]






![[Inbox 8.pdf]]


$$M_{11} = 1000+600C\sin(2\pi x_1)\sin(2\pi x_2)$$
$$M_{22} = 1000-600C\sin(2\pi x_1)\sin(2\pi x_2)$$
For the case $C=1.65$, so the $\max{M}/\min{M}\approx 20000$

![[New_C1.65.png]]

Linearization 1 (Previous work):
$$\frac{\partial x_i}{\partial s_j}\frac{\partial x_i}{\partial s_j} \to \frac{\partial x_i}{\partial s_j}\frac{(\partial x_i + \Delta x_i)}{\partial s_j} = \left(\frac{\partial x_i}{\partial s_j} \right)^2 + \frac{\partial  \Delta x_i}{\partial s_j}$$
Linearization 2:
$$\frac{\partial x_i}{\partial s_j}\frac{\partial x_i}{\partial s_j} \to \frac{(\partial x_i + \Delta x_i)}{\partial s_j}\frac{(\partial x_i + \Delta x_i)}{\partial s_j} = \left(\frac{\partial x_i}{\partial s_j} \right)^2 + 2\frac{\partial  \Delta x_i}{\partial s_j}$$
No Jacobian on $M$ so far

With C=1
![[C=1.png]]

With C=10
![[C=10.png]]


![[Inbox 7.pdf]]

Upwind is used for those first order derivative discretizations.


![[omega=0.1,C=1,Nx1=Nx2=50.png]]
![[omega=0.1,C=2,Nx1=Nx2=50.png]]
![[omega=0.1,C=10,Nx1=Nx2=50.png]]



Ok, last week, I thought that we can use Newton's method to solve the system instead. Now, what I found was the method mentioned in the paper can be though as Newton's method based on linear approximation of the original equation. Initially, we have:

$$F(x) = \sum_{\alpha=1}^3 A_{kj,\alpha}(x) D_{\alpha \alpha} x_k + R(x)=0$$
Then the Jacobian is:
$$J(x)_{km} = \frac{\delta F_k}{\delta x_m} = \sum_{\alpha=1}^3 \frac{\partial A_{kj,\alpha}(x)}{\partial x_m} D_{\alpha \alpha} x_k + A_{kj,\alpha}(x) D_{\alpha \alpha} \delta_{km} + \frac{\partial R(x)_k}{\partial x_m}$$
If we are assuming $\frac{\partial A_{kj,\alpha}(x)}{\partial x_m}$ and $\frac{\partial R(x)_k}{\partial x_m}$ is small (which is not true for case 1, and is exactly 0 for other cases), then the Jacobian is:
$$J(x)_{km} \approx \frac{\delta F_k}{\delta x_m} = \sum_{\alpha=1}^3 A_{kj,\alpha}(x) D_{\alpha \alpha} \delta_{km} $$
Then Newton's method gives:
$$x^{n+1} = x^n - \frac{F(x^n)}{J(x^n)}$$
Or:
$$J(x^n)\Delta x^n = -F(x^n)$$
Which turns out to be our iterations:
$$\sum_{\alpha=1}^3 A_{kj,\alpha}(x^n) D_{\alpha \alpha} \Delta x^n = -\text{Res}(x^n)$$
I also tried to use Damped Newton's method. It stops after some iterations because we already achieved minimum residual along the line search without residual converging to zero. So, I think the ultimate criminal is that the approximation is indeed not accurate enough: we need to introduce additional variational terms.


## Mar 14 Update:

Idea:
Iterative method bicgstabl with ilu preconditioner: 100\*100 grid --> 20,000\*20,000 matrix system converges in a couple of iterations. 

If we can get rid of x dependency of matrix A, then we can do just ilu once and solve all the iteration almost for free. 

Without preconditioner: Takes more than 10 sec 5e-4 residual 
With preconditioner: Takes about 0.4 sec to achieve 1e-6 residual



## Mar 7 Update: Better iteration
$$x^{n+1}=x^n+\omega \delta x^{**},\;\;\;\;\;\; \begin{cases}\mathcal{A}(\tilde{x}^{n+1}) (\tilde{x}^{n+1}+\delta x^{**})=\mathcal{R}(\tilde{x}^{n+1}) & \text{somehow works better}\\
\mathcal{A}(\tilde{x}^{n+1}) (x^n+\delta x^{**})=\mathcal{R}(\tilde{x}^{n+1})
\end{cases}$$
Where $$\tilde{x}^{n+1}=x^n+\omega \delta x^*,\;\;\;\;\;\;\mathcal{A}(x^n)(x^n+\delta x^*) = \mathcal{R}(x^n)$$

Make metric tensor more variant in magnitude:
$$M_{ii} = 40000(1+15Cx_i)^{-2}$$
![[convergence_comparison_1.png]]
![[convergence_comparison_2.png]]
![[convergence_comparison_3.png]]

## Feb 13 Update: Update Correction on Curved Boundary
For a curved boundary, allowing the boundary points to slide on the tangent plane— a local approximation of the actual boundary—can introduce truncation errors. For smooth boundaries, this error is typically small. However, for boundaries with high curvature, the deviation of the boundary points from the actual boundary can be significant.

Therefore, after each iteration, a correction process is required to ensure that updated boundary points remain on the predefined actual boundary. Since the actual boundary is typically defined numerically as a discrete set of points, an implicit linear interpolation scheme for parametrized curves is employed during correction.  

To correct an updated boundary point $p^*=(x^*,y^*)$ on the boundary, a reasonable choice would be the shortest projection of the updated boundary point on all of the interpolated segments. Let's look at an arbitrary segment with end points $p_i=(x_i,y_i)$ and $p_{i+1}=(x_{i+1},y_{i+1})$ from the array of actually boundary. Then we can have two vectors $u_i=p^*-p_i$ and $v_i=p_{i+1}-p_i$. Then the projection factor would be $$r=\frac{u_i \cdot v_i}{\Vert v_i\Vert^2}$$Notice that if $r<0$, then the orthogonal projection is before $p_i$. Meanwhile if $r>1$, then the orthogonal projection is after $p_{i+1}$. Therefore, we have to use the endpoints as the closest point on the segment by calling $r:=\max(0,\min(1,r))$.

Then the location of the projection on the segment $i$ would be $p^*_i=(x_i^*,y_i^*)=p_i+r*v_i$. Then the distance would be $\Vert p_i^* - p^*\Vert$. Therefore, we can compare all the segments on the boundary and decide which is closest segment which is the closest projection point.

## Feb 7 Update: Curved boundary
For curved boundaries, the Dirichlet boundary condition at left boundary becomes:
$$x_1^{n+1}(1,:)n_1 + x_2^{n+1}(1,:)n_1=x_1^{n}(1,:)n_1 + x_2^{n}(1,:)n_1$$
The local normal vector can be computed from the local tangent vector. The local tangent vector are approximated from discrete data at the previous iteration $t(1,j)=x^n(1,j+1)-x^n(1,j-1)$ with normalization.
This will results a 2-points stencil with right hand sides depending on the previous iteration.

For Neumann boundary condition, we use one-sided second order accurate derivative approximation, $\frac{d}{ds_\gamma}(x(1,:)\cdot t)\approx\frac{1}{2\sigma}[3(x(1,:)\cdot t)-4(x(2,:)\cdot t)+x(3,:)\cdot t)]=0$ which leads to $(x(1,:)\cdot t)=\frac{1}{3}[4(x(2,:)\cdot t)-(x(3,:)\cdot t)]$ , or:
$$x_1(1,:)t_1+x_2(1,:)t_2 - \frac{1}{3}\left[4(x_1(2,:)t_1+x_2(2,:)t_2)- (x_1(3,:)t_1+x_2(3,:)t_2) \right]=0$$
This will results a 6-points stencil with 0 on the right hand side. Alternative, we can use the interior information from previous iteration and put them on the right hand side.


#### Iteration animation (with 40 points in each direction, relaxation $\omega=0.2$)

![[example13.gif]]

To-do list:
- [ ] Implement non-smooth curved boundary
- [ ] Implement solvers with approximate cost function
- [ ] Improve code
	- [ ] Modulization
	- [ ] Create config - calculation - output isolation
	- [ ] Documentation
- [ ] Result analysis
	- [ ] Residual analysis
	- [ ] Relaxation vs convergence


---

![[Notes.pdf]]


#### Improvement
To resolve the Neumann boundary conditions $\frac{\partial x_2}{\partial s_1}=0$ on the left and right boundaries, and $\frac{\partial x_1}{\partial s_2}=0$ on the top and bottom boundaries. We can setup ghost points:

$x_2|_{s_1=-\sigma_1} := x_2|_{s_1=\sigma_1}$ on the left boundary
$x_2|_{s_1=1+\sigma_1} := x_2|_{s_1=1-\sigma_1}$ on the right boundary
$x_1|_{s_2=-\sigma_2} := x_1|_{s_2=\sigma_2}$ on the bottom boundary
$x_1|_{s_2=1+\sigma_2} := x_1|_{s_2=1-\sigma_2}$ on the top boundary

So that the central difference at $s_1=0$, $s_1=1$, $s_2=0$, $s_2=1$, will be zero respectively. This will resolve the Neumann boundary condition with second order accuracy.

To impose those relations in the matrix system, we will modify the approximated second order derivatives in those rows associated with boundary points. For example, $\partial^2 x_2/\partial s_1^2$ on the left boundary would become:
$$\frac{\partial^2 x_2}{\partial s_1^2}\approx \frac{x_2|_{s_1=-\sigma_1} - 2x_2|_{s_1=0} + x_2|_{s_1=\sigma_1}}{\sigma_1^2} = 2\frac{x_2|_{s_1=\sigma_1} - x_2|_{s_1=0}}{\sigma_1^2}$$
This is more accurate that directly enforcing the boundary term $x_2|_{s_1=0} = x_2|_{s_1=\sigma_1}$. We will see that the results will be smoother even without orthogonalization.


#### Implementation
~~~
For iter from 1 to max_ter:
	[A, b, residual] = AssembleMatrixSystem(x)
	x_new = A \ b
	dx = x_new - x
	x = x + omega*dx
	Break if residual < tol
End
~~~
~~~
AssembleMatrixSystem:
	Get points id
	Compute coefficient associated with each point
	For each component of s (s1, s2):
		For each geometrix position (top, bottom, left, right, interior)
			A(associated_row, associated_row + stencil_shift) += coefficient * difference_coefficient;
		end
	end
	Compute b given by coefficient
	Let A's rows constrained by Dirichlet condition be identity
	Let b's rows constrained by Dirichlet condition be 0 or 1
	Compute residual
~~~


#### Assembled matrices (5x5 points mesh)

Alternative
![[ExampleAlterMatrix.png]]

Approximated
![[ExampleApproxMatrix.png]]


#### Iteration animation (with 40 points in each direction, relaxation $\omega=0.2$)

![[example11.gif]]
![[example21.gif]]

![[example12.gif]]
![[example22.gif]]



Questions:
1. The Naumann condition for curved boundary commits variational crime: boundaries will deform as iteration goes on
	1. Is it reasonable to directly enforce boundary points onto the boundary curves after each iteration?
3. Relevant papers about metric conforming structured mesh, multi-block, etc.