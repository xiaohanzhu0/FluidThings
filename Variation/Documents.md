
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