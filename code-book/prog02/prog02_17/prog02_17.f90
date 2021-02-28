!##############################################################################
! PROGRAM polynomials
!
! ## Polynomial interpolation using different kinds of nodes
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
program polynomials

    use toolbox

    implicit none
    integer, parameter :: n = 10, nplot = 1000
    real*8 :: xi(0:n), yi(0:n)
    real*8 :: xplot(0:nplot), yplot(0:nplot), yreal(0:nplot)

    ! get equidistant plot nodes and Runge's function
    call grid_Cons_Equi(xplot, -1d0, 1d0)
    yreal = 1d0/(1d0+25*xplot**2)
    call plot(xplot, yreal, legend='Original')

    ! equidistant polynomial interpolation
    call grid_Cons_Equi(xi, -1d0, 1d0)
    yi = 1d0/(1d0+25*xi**2)
    yplot = poly_interpol(xplot, xi, yi)
    call plot(xplot, yplot, legend='Equidistant')

    ! Chebyshev polynomial interpolation
    call grid_Cons_Cheb(xi, -1d0, 1d0)
    yi = 1d0/(1d0+25*xi**2)
    yplot = poly_interpol(xplot, xi, yi)
    call plot(xplot, yplot, legend='Chebyshev')

    ! execute plot
    call execplot()

end program
