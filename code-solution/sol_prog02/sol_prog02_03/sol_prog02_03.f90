!##############################################################################
! PROGRAM GoldenSearch
!
! ## Compute the market equilibrium using spline interpolation
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog02_03m.f90"

program GoldenSearch

    use toolbox
    use globals

    implicit none

    ! declaration of variables
    real*8 :: min_global, fmin_global
    integer :: i_global, i

    ! set up n subintervals over [0,5]
    call grid_Cons_Equi(x, x_l, x_u)

    ! use the function minimize find minimum value in each interval
    do i = 1, n
        minimum_x(i) = minimize(x(i-1), x(i))
    enddo

    ! calculate function values
    fmin = minimum_x*cos(minimum_x**2d0)

    ! locate global minimum
    i_global = minloc(fmin, 1)
    min_global = minimum_x(i_global)
    fmin_global = min_global*cos(min_global**2d0)

    ! print output to the screen
    write(*,'(a)')'The global minimum is lcated at'
    write(*,'(2(a,f10.6))')'x = ', min_global, '    y = ', fmin_global

    ! set up data for plot
    call grid_Cons_Equi(xplot, x_l, x_u)
    yplot = xplot*cos(xplot**2d0)

    ! initialize plot
    call plot(xplot, yplot)

    ! execute plot
    call execplot(xlabel='x', ylabel='x cos(x^2)')

end program
