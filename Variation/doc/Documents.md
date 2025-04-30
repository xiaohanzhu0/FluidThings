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