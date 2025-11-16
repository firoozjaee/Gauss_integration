program test_gauss_mod
    use Gauss_integration
    implicit none
    real(kind=8) :: result
    real(kind=8) :: a(2) = [-1.0d0, -1.0d0]
    real(kind=8) :: b(2) = [1.0d0 ,  1.0d0]
    integer :: ip
    
    ! Set user's custom function
    call set_integrand_function(my_custom_function)
    
    ! Simple usage
    call simple_integrate(2, 5, a, b, result)
    print '(A,F20.8)', "2D Gaussian integral:", result
    
    ! Even simpler usage with uniform limits
    do ip = 100, 100, 100
        call integrate_uniform(2, ip, -1.0d0, 1.0d0, result)
        print '(A,F20.8,I6)', "2D Uniform integral :", result, ip
    end do

contains

    ! Definition of user's custom function
    real(kind=8) function my_custom_function(x, n)
        integer, intent(in) :: n
        real(kind=8), intent(in) :: x(n)
        ! Here you can define any custom function
        ! my_custom_function = 1.0d0  ! Constant function 1
        ! my_custom_function = x(1)   ! Linear function
        
        my_custom_function = exp(-sum(x**2))  ! Gaussian function
        ! my_custom_function = x(1)**2 + x(2)**2  ! Quadratic function
        ! my_custom_function = sin(x(1)**2) * cos(x(2)**2)  ! Trigonometric function
    end function my_custom_function

end program test_gauss_mod