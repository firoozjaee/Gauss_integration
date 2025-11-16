module Gauss_integration
!-----------------------------------------------------
!
!   This code is produced by Ali Rahmani Firoozjaee
!   Civil Engineering Department
!   Noushirvani University of Technology
!   
!-----------------------------------------------------
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    
    ! Only these three subroutines are accessible to the user
    public :: simple_integrate, integrate_uniform, set_integrand_function
    
    ! All other subroutines remain hidden
    private
    private :: nD_gauss_quadrature, gauss_legendre, legendre_polynomial
    private :: sort_points_weights, integrand_nd, default_integrand

    ! Interface definition for user's integrand function
    interface
        function integrand_interface(x, n) result(f)
            import dp
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp) :: f
        end function
    end interface    
    
    ! Pointer to user-defined integrand function
    procedure(integrand_interface), pointer :: user_integrand
    
contains

    ! Subroutine to set user's custom integrand function
    subroutine set_integrand_function(func)
        procedure(integrand_interface) :: func
        user_integrand => func
    end subroutine

    ! Default integrand function if user doesn't set one
    real(dp) function default_integrand(x, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        default_integrand = 1.0_dp
    end function default_integrand

    ! Simple interface for user - only essential parameters
    subroutine simple_integrate(n_vars, n_points, a, b, result)
        integer, intent(in) :: n_vars, n_points
        real(dp), intent(in) :: a(n_vars), b(n_vars)
        real(dp), intent(out) :: result
        real(dp), allocatable :: points(:), weights_gauss(:), x_transformed(:)
        real(dp) :: initial_wei, initial_detj
        character(len=100) :: filename
        integer :: i, ierr
        
        ! Memory allocation
        allocate(points(n_points), weights_gauss(n_points), x_transformed(n_vars))
        
        ! Generate or read Gauss points
        write(filename, '("gauss_points_", I0, ".txt")') n_points
        open(unit=10, file=filename, status='old', action='read', iostat=ierr)
        if (ierr == 0) then
            read(10, *) ! skip n
            do i = 1, n_points
                read(10, *) points(i), weights_gauss(i)
            end do
            close(10)
            print *, "Loaded Gauss points from file: ", filename
        else
            call gauss_legendre(n_points, points, weights_gauss)
            ! print *, "Generated Gauss points for n =", n_points
        end if
        
        ! Initialization and main call
        result = 0.0_dp
        x_transformed = 0.0_dp    
        initial_wei = 1.0_dp
        initial_detj = 1.0_dp
        
        call nD_gauss_quadrature(n_vars, n_vars, n_points, a, b, points, weights_gauss, &
                                x_transformed, initial_wei, initial_detj, result)
        
        ! Memory deallocation
        deallocate(points, weights_gauss, x_transformed)
    end subroutine simple_integrate

    ! Even simpler interface for uniform integration limits
    subroutine integrate_uniform(n_vars, n_points, a, b, result)
        integer, intent(in) :: n_vars, n_points
        real(dp), intent(in) :: a, b  ! Same limits for all dimensions
        real(dp), intent(out) :: result
        real(dp) :: a_vec(n_vars), b_vec(n_vars)
        
        a_vec = a
        b_vec = b
        call simple_integrate(n_vars, n_points, a_vec, b_vec, result)
    end subroutine integrate_uniform

    ! =========================================================================
    ! Internal subroutines - user should not call these directly
    ! =========================================================================

    ! Main integrand function that uses either user function or default
    real(dp) function integrand_nd(x, n)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        if (associated(user_integrand)) then
            integrand_nd = user_integrand(x, n)
        else
            integrand_nd = default_integrand(x, n)
        end if
    end function integrand_nd

    ! Recursive N-dimensional Gauss quadrature subroutine
    recursive subroutine nD_gauss_quadrature(ndim, n_vars, n_points, a, b, points, weights_gauss, &
                                            x_transformed, current_wei, current_detj, integral)
        integer, intent(in) :: ndim, n_vars, n_points
        real(dp), intent(in) :: a(ndim), b(ndim)
        real(dp), intent(in) :: points(n_points), weights_gauss(n_points)
        real(dp), intent(inout) :: x_transformed(ndim)
        real(dp), intent(inout) :: current_wei, current_detj
        real(dp), intent(inout) :: integral
        real(dp) :: new_wei, new_detj, temp_integral
        integer :: i, current_dim
        
        current_dim = ndim - n_vars + 1
        
        if (n_vars == 1) then
            ! Base case: last dimension
            do i = 1, n_points
                x_transformed(current_dim) = ((b(current_dim) - a(current_dim)) * points(i) + &
                                            (a(current_dim) + b(current_dim))) / 2.0_dp
                new_wei = weights_gauss(i)
                new_detj = (b(current_dim) - a(current_dim)) / 2.0_dp
                
                integral = integral + (current_wei * new_wei) * (current_detj * new_detj) * &
                                   integrand_nd(x_transformed, ndim)
            end do
        else
            ! Recursive case: remaining dimensions
            do i = 1, n_points
                x_transformed(current_dim) = ((b(current_dim) - a(current_dim)) * points(i) + &
                                            (a(current_dim) + b(current_dim))) / 2.0_dp
                new_wei = weights_gauss(i)
                new_detj = (b(current_dim) - a(current_dim)) / 2.0_dp
                
                new_wei = current_wei * new_wei
                new_detj = current_detj * new_detj

                temp_integral = 0.0_dp
                call nD_gauss_quadrature(ndim, n_vars-1, n_points, a, b, &
                                       points, weights_gauss, x_transformed, &
                                       new_wei, new_detj, temp_integral)
                integral = integral + temp_integral
            end do
        end if
    end subroutine nD_gauss_quadrature

    ! Generate Gauss-Legendre points and weights
    subroutine gauss_legendre(n, points, weights)
        integer, intent(in) :: n
        real(dp), intent(out) :: points(n), weights(n)
        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
        real(dp) :: x, f, dx, pn(0:n), pn_deriv
        integer :: i, iterations
        real(dp) :: tolerance = 1.0e-15_dp
        
        do i = 1, n
            x = cos(pi * (i - 0.25_dp) / (n + 0.5_dp))
            iterations = 0
            
            do while (iterations < 100)
                call legendre_polynomial(n, x, pn)
                f = pn(n)
                
                if (abs(x) == 1.0_dp) then
                    pn_deriv = n * (n + 1) * 0.5_dp * x**(n+1)
                else
                    pn_deriv = n * (x * pn(n) - pn(n-1)) / (x**2 - 1.0_dp)
                end if
                
                dx = f / pn_deriv
                x = x - dx
                iterations = iterations + 1
                if (abs(dx) < tolerance) exit
            end do
            
            points(i) = x
            weights(i) = 2.0_dp / ((1.0_dp - x**2) * pn_deriv**2)
        end do
        
        call sort_points_weights(n, points, weights)
    end subroutine gauss_legendre

    ! Compute Legendre polynomial values
    subroutine legendre_polynomial(n, x, p)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p(0:n)
        integer :: k
        
        p(0) = 1.0_dp
        if (n >= 1) p(1) = x
        
        do k = 2, n
            p(k) = ((2.0_dp * k - 1.0_dp) * x * p(k-1) - (k - 1.0_dp) * p(k-2)) / k
        end do
    end subroutine legendre_polynomial

    ! Sort points and weights in descending order
    subroutine sort_points_weights(n, points, weights)
        integer, intent(in) :: n
        real(dp), intent(inout) :: points(n), weights(n)
        real(dp) :: temp_p, temp_w
        integer :: i, j
        
        do i = 1, n-1
            do j = i+1, n
                if (points(i) < points(j)) then
                    temp_p = points(i)
                    points(i) = points(j)
                    points(j) = temp_p
                    
                    temp_w = weights(i)
                    weights(i) = weights(j)
                    weights(j) = temp_w
                end if
            end do
        end do
    end subroutine sort_points_weights

end module Gauss_integration