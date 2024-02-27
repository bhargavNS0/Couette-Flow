# Couette-Flow
Numerical Analysis of Couette flow using Crank-Nicolson scheme and Explicit scheme in Finite Difference Method.

Couette flow - It is the flow between two parallel plates, where one plate is held stationary and other is moving at constant velocity.
Mathematically, the velocity profile of Couette flow can be described by the linear equation: u(y) = (y/H)*U

The timestep for the Crank-Nickolson and Explicit scheme is taken based on CFL criteria i.e., dt/(Re*dy^2) <= 1/2
With the chosen values of Reynolds number and dy (deltaY), dt (deltaT) is calculated.

The tridiagonal matrix is solved using TDMA (Tri-Diagonal Matrix Algorithm) also known as Thomas Algorithm.


