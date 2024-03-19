# Finite Element Method for Solving Poisson's Equation

This project is a part of the course on Differential and Difference Equations. It focuses on solving Poisson's equation using the Finite Element Method (FEM).

## Introduction to the Problem: Gravitational Potential

The problem at hand involves solving for the gravitational potential Φ(x), which is governed by the following second-order differential equation:

$$\frac{d^2\Phi}{dx^2} = 4\pi G\rho(x)$$

In this equation:
- G represents the gravitational constant.
- ρ(x) denotes the density distribution of the mass.

We also have the following boundary conditions for Φ(x):
- At x = 0, Φ(0) = 5
- At x = 3, Φ(3) = 4

The density function ρ(x) is defined piecewise as follows:

$$
\rho(x) = 
\begin{cases} 
0 & \text{for } x \in [0, 1] \\
1 & \text{for } x \in (1, 2] \\
0 & \text{for } x \in (2, 3] 
\end{cases}
$$

## Resulting function

![graph (2)](https://github.com/stylub/Solving_Poisson_Equation_with_FEM/assets/47119994/31f0239d-aa58-42fe-bad0-12ae68c37350)
