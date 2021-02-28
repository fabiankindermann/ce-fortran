!##############################################################################
! PROGRAM CosFunc
!
! ## Interpolate Cosine Function linearly and with a spline
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
program CosFunc

    use toolbox

    implicit none

    ! set parameters
    integer, parameter :: n = 10
    integer, parameter :: n_plot = 1000
    real*8, parameter :: pi = 3.1415926535d0

    ! declaration of variables
    real*8 :: x(0:n), y(0:n), coeff(n+3)
    real*8 :: x_plot(0:n_plot), y_plot(3, 0:n_plot), y_real(0:n_plot)
    real*8 :: m(0:n), t(0:n), varphi(0:n_plot)
    integer :: i, j, il(0:n_plot), ir(0:n_plot)

    ! get nodes and data for plotting
    call grid_cons_Equi(x_plot, 0d0, 2d0*pi)
    y_real = cos(x_plot)

    ! get nodes and data for interpolation
    call grid_Cons_Equi(x, 0d0, 2d0*pi)
    y = cos(x)

    ! compute slope and intercept for each subinterval
    do i = 1, n

        m(i) = (y(i) - y(i-1))/(x(i) - x(i-1))
        t(i) = y(i) - m(i)*x(i)

    enddo

    ! piecewise linear interpolation by hand
    do j = 0, n_plot

        ! get the subinterval for a certain entry from xplot
        i = int(grid_Inv_Equi(x_plot(j), 0d0, 2d0*pi, n)) + 1

        ! k should be smaller than the maximum number of subintervals n
        i = min(i, n)

        ! interpolate function value
        y_plot(1, j) = m(i)*x_plot(j) + t(i)

    enddo

    ! piecewise linear interpolation with subroutine from toolbox
    call linint_Equi(x_plot, 0d0, 2d0*pi, n, il, ir, varphi)
    y_plot(2, :) = varphi*y(il) + (1d0-varphi)*y(ir)

    ! cubic spline interpolation
    call spline_interp(y, coeff)
    y_plot(3, :) = (/(spline_eval(x_plot(j), coeff, 0d0, 2d0*pi), j = 0, n_plot)/)

    ! print output
    write(*,'(a)') 'Approximation Error'
    write(*,'(a)') '-------------------'
    write(*,'(a, f15.10)') 'Piecewise Linear (Linear EQ)', maxval(y_plot(1, :)-y_real)
    write(*,'(a, f15.10)') 'Piecewise Linear (Toolbox)', maxval(y_plot(2, :)-y_real)
    write(*,'(a, f15.10)') 'Cubic Spline', maxval(y_plot(3, :)-y_real)

    ! initialize plot
    call plot(x_plot, y_plot(1, :), legend='Piecewise Linear (Linear EQ)')
    call plot(x_plot, y_plot(2, :), legend='Piecewise Linear (Toolbox)')
    call plot(x_plot, y_plot(3, :), legend='Cubic Spline')
    call plot(x_plot, y_real, legend='Cos(x)')

    ! execute plot
    call execplot(xlabel='x', ylabel='Cos(x)')

end program
