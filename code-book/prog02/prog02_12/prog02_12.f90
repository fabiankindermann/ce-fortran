!##############################################################################
! PROGRAM NewtonCotes
!
! ## Numerical integration with Newton-Cotes formulas
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
program NewtonCotes

    implicit none
    integer, parameter :: n = 10
    real*8, parameter :: a = 0d0, b = 2d0
    real*8 :: h, x(0:n), w(0:n), f(0:n)
    integer :: i

    ! calculate quadrature nodes
    h = (b-a)/dble(n)
    x = (/(a + dble(i)*h, i=0,n)/)

    ! get weights
    w(0) = h/2d0
    w(n) = h/2d0
    w(1:n-1) = h

    ! calculate function values at nodes
    f = cos(x)

    ! Output numerical and analytical solution
    write(*,'(a,f10.6)')' Numerical:  ',sum(w*f, 1)
    write(*,'(a,f10.6)')' Analytical: ',sin(2d0)-sin(0d0)

end program
