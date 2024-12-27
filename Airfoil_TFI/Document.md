#### Simplest Case: Unidirectional Interpolation ####

1. Interpolation is performed only in the $\xi^1 \to x^1$ direction.  
2. Zero-th order mapping of data points is applied, i.e., $P_l^i = 0$.  
3. Only boundary data points are used for interpolation, with $L^i = 2$.  
4. Linear Lagrange interpolation is employed for the blending term, with $\alpha_1 = 1 - \xi$ and $\alpha_2 = \xi$.  

Under these conditions, the mapping in Equation (5.1) simplifies to:

$$
\mathbf{P}_2(\boldsymbol{\xi}) = (1 - \xi^1) \mathbf{A}^1(0, \xi^2) + \xi^1 \mathbf{A}^1(1, \xi^2),
$$

where $\mathbf{A}^i(\xi^1, \xi^2)$ represents the map $\xi \to x$.  

In this case:
- $\xi^1 = 0$ or $\xi^1 = 1$, and  
- $\mathbf{A}^1(\xi^1, \xi^2)$ is defined only at a finite set of $\xi^2$ values.  

Thus:  
- Only the boundary data points at the left and right boundaries are used for interpolation.  
- The boundary data points at the top and bottom boundaries are not utilized.  
- For simplicity, we can omit the superscript. If the i-th variable is a constant, then we know it has a superscipt of i.

As a result, the interpolated value varies linearly in the $\xi^1$ direction.

For implementation, the data at the y boundary is given. $\mathbf{A}^1(0,\xi^2)$ is equivalent to `jneg`. If we want to draw a set of nodes equi-spaced in $\xi^1$, and with $\xi^2$ inside `jneg` and `jpos` . Then we have
```s = (1-xi1).*reshape(jneg,1,Nj,2) + xi1.*reshape(jpos,1,Nj,2);```
s will be a 3D array with Ni\*Nj\*Nd

The generation result is shown below:
![[UniLinear.png]]

#### Follow on: Bidirectional interpolation
Every condition is the same as the unidirectional case, except the first one: now interpolation in both $\xi^1 \to x^1$ and $\xi^2 \to x^2$ direction. Under these conditions, the general formula given in Equation (5.19) simplifies to:
$$
\begin{align}
\mathbf{P}(\boldsymbol{\xi}) =& (1-\xi^2)\mathbf{A}(\xi^1,0) + \xi^2\mathbf{A}(\xi^1,1)+(1-\xi^1)\mathbf{A}(0,\xi^2) + \xi^1\mathbf{A}(1,\xi^2) \\
&- \left[ (1-\xi^1)(1-\xi^2)\mathbf{A}(0,0) + \xi^1\xi^2 \mathbf{A}(1,1) + \xi^1(1-\xi^2)\mathbf{A}(1,0) + (1-\xi^1)\xi^2 \mathbf{A}(0,1) \right]
\end{align}
$$
Which is consistent with the formula given in the Wikipedia page.

![[demonstration.pdf]]

The generation result is shown below:
![[BiLinear.png]]

As observed, the results from the unidirectional and bidirectional methods are quite similar. In fact, the generated interior cells are identical to machine precision. This is expected since the boundary data points in the $i$-direction are evenly spaced. To highlight the differences, let us modify the boundary data points in the $i$-direction:

Curved Boundary:
![[BiLinear_curved.png]]

Non equi-spaced:
![[BiLinear_nonequispaced.png]]

Curved and non equi-spaced:
![[BiLinear_curved_nonequispaced.png]]


#### Address on the Orthogonality####
The orthogonality of interior cells depends heavily on the boundary definitions. If the boundary data points are equi-spaced, the resulting interior grid will naturally exhibit orthogonality. 

My guess: However, if the boundary data points are not equi-spaced, it seems like that there is no feasible way to keep orthogonality. More advanced methods are need, such many of the optimization-based method involving solving differential equations.



