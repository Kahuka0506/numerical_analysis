# Numerical analysis

Numerical analysis library

## Linear algebra  
[Linear algebra](./linear_algebla/README.md)
- Direct method `direct_method.hpp`
    - Gauss Elimination, Matrix Rank `gauss_elimination.cpp` 
    - LU decomposition, Inverse matrix, Linear equations `LU.cpp` 
    - Cholesky

- Iterative method `iterative_method.hpp`  
    - Stationary method `linear_equation.cpp` Jacobi, Gauss-Seidel, SOR    
    - Nonstationary method `CG.cpp` CG, BCG, CGS, BISGSTAB  


- Eigen value `eigen.hpp`
    - Power method, Power method shift, Inverse power method `power_method.cpp`
    - QR, Howseholder




## Differential Equations
[Differential Equations](./differential_equations/README.md)
- Ordinary differential equations(ODE)
    - Runge kutta
- Partial differential equations(PDE)
    - Elliptic (B^2-AC < 0)  Poisson equation    
    - Parabolic (B^2-AC = 0)  Diffusion equation  
    Explicit,Implicit,Crank-Nicolson  
    - Hyperbolic (B^2 - AC > 0)   Wave equation    


## Interpolation
[Interpolation](./interpolation/README.md)
- Lagrange interpolation `lagrange_interpolation.cpp` 
- Spline interpolation `spline_interpolation.cpp` 
- Least squares `least_squares.cpp`




## Integration
[Integration](./integration/README.md)
- Middle method
- Trapezoid method
- Simpson method
- Simpson_3_8 method



