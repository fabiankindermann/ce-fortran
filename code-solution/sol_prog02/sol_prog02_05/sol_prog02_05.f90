!##############################################################################
! PROGRAM Integral
!
! ## Compute the integal for of exp(-x) and |x|^0.5 numerically
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
! #VC# VERSION: 1.0  (16 December 2019)
!
!##############################################################################
program integral

    use toolbox

    implicit none

    ! number of quadrature nodes
    integer, parameter :: n = 10

    ! switch whether to use symmetry to compute the integral for |x|^0.5
    logical :: symmetry = .false.

    ! variables for nodes, weights and function values
    integer :: i
    real*8 :: h, a, b
    real*8 :: x_int(0:n), w_trapez(0:n), x_gauss(0:n), w_gauss(0:n)
    real*8 :: f_trapez_exp(0:n), f_gauss_exp(0:n), f_trapez_sqrtabs(0:n), f_gauss_sqrtabs(0:n)
    real*8 :: analytical_exp, trapez_exp, simpson_exp, gauss_exp
    real*8 :: analytical_sqrtabs, trapez_sqrtabs, simpson_sqrtabs, gauss_sqrtabs

    ! set lower and upper bound for interval
    a = -1d0
    b = 1d0

    ! get equidistant nodes
    call grid_Cons_Equi(x_int, a, b)

    ! compute interval length
    h = (b-a)/dble(n)

    ! analyticalal solution
    analytical_exp = exp(1d0) - exp(-1d0)

    ! compute weights for trapezoid rule
    w_trapez(0) = h/2d0
    w_trapez(n) = h/2d0
    w_trapez(1:n-1) = h

    ! use trapezoid rule for integral of exp(-x)
    f_trapez_exp = exp(-x_int)
    trapez_exp = sum(w_trapez*f_trapez_exp,1)

    ! use simpson rule for integral
    simpson_exp = exp(-a) + exp(-b)

    do i = 1, n-1
        simpson_exp = simpson_exp + 2d0*exp(-x_int(i))
    enddo

    do i = 0, n-1
        simpson_exp = simpson_exp + 4d0*exp(-(x_int(i)+x_int(i+1))/2d0)
    enddo

    simpson_exp = simpson_exp*h/6d0

    ! gauss-legendre nodes and weights
    call legendre(a, b, x_gauss, w_gauss)

    ! gauss-legendre for integral of exp(-x)
    f_gauss_exp = exp(-x_gauss)
    gauss_exp = sum(w_gauss*f_gauss_exp, 1)

    ! print output
    write(*,'(a, f15.9)') 'Analytical Solution (exp(-x)):', analytical_exp
    write(*,'(/a)') '                   Trapez             Simpson      Gauss-Legendre'
    write(*,'(a, 3f20.9)') 'Delta', trapez_exp - analytical_exp, simpson_exp - analytical_exp, &
        gauss_exp - analytical_exp

    ! if symmetry property is used for the intergal of |x|^0.5,
    ! we set new interval borders and compute new weights
    if(symmetry)then

        ! set lower and upper bound for interval
        a = 0d0
        b = 1d0

        ! compute interval length
        h = (b-a)/dble(n)

        ! get equidistant nodes
        call grid_Cons_Equi(x_int, a, b)

        ! compute weights for trapezoid rule
        w_trapez(0) = h/2d0
        w_trapez(n) = h/2d0
        w_trapez(1:n-1) = h

        ! gauss-legendre nodes and weights
        call legendre(a, b, x_gauss, w_gauss)

    endif

    ! use trapezoid rule for integral of |x|^0.5
    f_trapez_sqrtabs = abs(x_int)**0.5d0
    trapez_sqrtabs = sum(w_trapez*f_trapez_sqrtabs, 1)

    ! use simpson rule for integral
    simpson_sqrtabs = abs(a)**0.5d0 + abs(b)**0.5d0

    do i = 1, n-1
        simpson_sqrtabs = simpson_sqrtabs + 2d0*abs(x_int(i))**0.5d0
    enddo

    do i = 0, n-1
        simpson_sqrtabs = simpson_sqrtabs + 4d0*abs((x_int(i)+x_int(i+1))/2d0)**0.5d0
    enddo

    simpson_sqrtabs = simpson_sqrtabs*h/6d0

    ! gauss-legendre for integral of |x|^0.5
    f_gauss_sqrtabs = abs(x_gauss)**0.5d0
    gauss_sqrtabs = sum(w_gauss*f_gauss_sqrtabs, 1)

    ! if symmetry is used multiply temporary results
    ! by two to get actual result
    if(symmetry)then
        trapez_sqrtabs = 2d0*trapez_sqrtabs
        simpson_sqrtabs = 2d0*simpson_sqrtabs
        gauss_sqrtabs = 2d0*gauss_sqrtabs
    endif

    ! analyticalal solution
    analytical_sqrtabs = 4d0/3d0

    ! print output
    write(*,'(/a, f15.9)') 'Analytical Solution (|x|^0.5):', analytical_sqrtabs
    write(*,'(/a)') '                   Trapez             Simpson      Gauss-Legendre'
    write(*,'(a, 3f20.9)') 'Delta', trapez_sqrtabs - analytical_sqrtabs, simpson_sqrtabs - analytical_sqrtabs, &
        gauss_sqrtabs - analytical_sqrtabs

end program
