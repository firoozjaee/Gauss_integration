# Gauss Integration Module

A Fortran module for multi-dimensional numerical integration using Gaussian quadrature with Legendre polynomials.

## Features

- **Multi-dimensional integration**: Supports N-dimensional integrals
- **Gauss-Legendre quadrature**: High accuracy with minimal function evaluations
- **User-friendly interface**: Simple and uniform integration interfaces
- **Custom integrand functions**: Easy integration of user-defined functions
- **Automatic point generation**: Generates Gauss points or loads from file
- **Memory efficient**: Automatic memory management

## Module Structure

### Public Subroutines

- `simple_integrate(n_vars, n_points, a, b, result)`
  - General integration with different limits for each dimension
  - `n_vars`: Number of dimensions
  - `n_points`: Gauss points per dimension
  - `a, b`: Integration limits arrays
  - `result`: Integration result

- `integrate_uniform(n_vars, n_points, a, b, result)`
  - Simplified interface with uniform limits for all dimensions
  - `a, b`: Single values applied to all dimensions

- `set_integrand_function(func)`
  - Sets user-defined integrand function
  - `func`: User's custom function following the integrand interface

### Integrand Function Interface

User-defined functions must follow this signature:
```fortran
real(kind=8) function my_function(x, n)
    integer, intent(in) :: n
    real(kind=8), intent(in) :: x(n)
    ! Your function definition here
end function

! Compilation
# Compile both files together
gfortran -o program.exe main_program.f90 Gauss_integration.f90

# Or compile separately
gfortran -c Gauss_integration.f90
gfortran -c main_program.f90
gfortran -o program.exe main_program.o Gauss_integration.o