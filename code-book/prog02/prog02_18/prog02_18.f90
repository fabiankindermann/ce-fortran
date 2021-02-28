!##############################################################################
! PROGRAM piecewise_pol
!
! ## Piecewise polynomial interpolation (linear and cubic splines)
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
program piecewise_pol

    use toolbox

    implicit none
    integer, parameter :: n = 10, nplot = 1000
    real*8 :: xi(0:n), yi(0:n)
    real*8 :: xplot(0:nplot), yplot(0:nplot), yreal(0:nplot)
    real*8 :: varphi, coeff(n+3)
    integer :: ix, il, ir

    ! get nodes and data for interpolation
    call grid_Cons_Equi(xi, -1d0, 1d0)
    yi = 1d0/(1d0+25*xi**2)

    ! get nodes and data for plotting
    call grid_Cons_Equi(xplot, -1d0, 1d0)
    yreal = 1d0/(1d0+25*xplot**2)
    call plot(xplot, yreal, legend='Original')

    ! piecewise linear interpolation
    do ix = 0, nplot
        call linint_Equi(xplot(ix), -1d0, 1d0, n, il, ir, varphi)
        yplot(ix) = varphi*yi(il) + (1d0-varphi)*yi(ir)
    enddo
    call plot(xplot, yplot, legend='Piecewise linear')

    ! cubic spline interpolation
    call spline_interp(yi, coeff)
    do ix = 0, nplot
        yplot(ix) = spline_eval(xplot(ix), coeff, -1d0, 1d0)
    enddo
    call plot(xplot, yplot, legend='Cubic spline')

    ! execute plot
    call execplot(xlim=(/-1d0,1d0/))

end program
