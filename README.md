# Solving Poisson Equation with FEM
### Project made as a part of course Differential and Difference Equations
## Gravitational Potential

The second-order differential equation for the gravitational potential Φ(x) is:

$$\frac{d^2\Phi}{dx^2} = 4\pi G\rho(x)$$

where:
- G is the gravitational constant,
- ρ(x) is the density of the mass distribution.

Boundary conditions for Φ(x) are as follows:
- Φ(0) = 5
- Φ(3) = 4

The density function ρ(x) is defined as:

$$
\rho(x) = 
\begin{cases} 
0 & \text{for } x \in [0, 1] \\
1 & \text{for } x \in (1, 2] \\
0 & \text{for } x \in (2, 3] 
\end{cases}
$$


![graph (2)](https://github.com/stylub/Solving_Poisson_Equation_with_FEM/assets/47119994/31f0239d-aa58-42fe-bad0-12ae68c37350)
