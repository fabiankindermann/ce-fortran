!##############################################################################
! PROGRAM GaussLegendre
!
! ## Numerical integration with Gauss-Legendre quadrature
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
! #VC# VERSION: 1.0  (23 January 2018)
!
!##############################################################################
program GaussLegendre

    use toolbox

    implicit none
    integer, parameter :: n = 10
    real*8, parameter :: a = 0d0, b = 2d0
    real*8 :: x(0:n), w(0:n), f(0:n)

    ! calculate nodes and weights
    call legendre(0d0, 2d0, x, w)

    ! calculate function values at nodes
    f = cos(x)

    ! Output numerical and analytical solution
    write(*,'(a,f10.6)')' Numerical:  ',sum(w*f, 1)
    write(*,'(a,f10.6)')' Analytical: ',sin(2d0)-sin(0d0)

end program
