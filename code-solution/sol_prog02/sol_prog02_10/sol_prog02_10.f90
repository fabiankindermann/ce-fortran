!##############################################################################
! PROGRAM Monopoly
!
! ## Interpolate Cosine Function linearly and with a spline
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog02_10m.f90"

program Monopoly

    use globals
    use toolbox

    implicit none

    ! declaration of variables
    logical :: check
    integer :: ir, ia, ix(2)
    real*8 :: x_in(2), fret

    ! initialize profit matrix
    G(0, :) = (/11.5d0, 70.9d0, 98.3d0, 93.7d0/)
    G(1, :) = (/31.1d0, 82.5d0, 101.9d0, 89.3d0/)
    G(2, :) = (/18.7d0, 62.1d0, 73.5d0, 52.9d0/)
    G(3, :) = (/-25.7d0, 9.7d0, 13.1d0, -15.5d0/)

    ! initialize coefficients for interpolation to zero
    coeff_G = 0d0

    ! set up grid for interpolation
    call grid_Cons_Equi(pa, p_min, p_max)
    call grid_Cons_Equi(pr, p_min, p_max)

    ! set up grid for evaluation of the profit function
    call grid_Cons_Equi(pa_plot, p_min, p_max)
    call grid_Cons_Equi(pr_plot, p_min, p_max)

    ! interpolate coefficients
    call spline_interp(G, coeff_G)

    ! two-dimensional interpolation for profit functions
    do ir = 0, Nplot
        do ia = 0, Nplot
            G_plot(ir, ia) = spline_eval((/pr_plot(ir), pa_plot(ia)/), coeff_G, (/p_min, p_min/), (/p_max, p_max/))
        enddo
    enddo

    ! plot the profit function using spline interpolation
    call plot3d(pr_plot, pa_plot, G_plot, xlabel = 'p_R', ylabel = 'p_A', zlabel = 'G(p_R, p_A)')

    ! locate maximum of the interpolated profit function
    ix = maxloc(G_plot)-1

    ! print result to the screen
    write(*,'(a)') '                       P_R         P_A      Profit'
    write(*,'(a)') '--------------------------------------------------'
    write(*,'(a, 3f12.6)') 'Interpolation:', pr_plot(ix(1)), pa_plot(ix(2)), G_plot(ix(1), ix(2))

    ! initialize starting value for fminsearch
    x_in = 0.5d0

    call fminsearch(x_in, fret, (/p_min, p_min/), (/p_max, p_max/), profit)

    ! print result to the screen
    write(*,'(a, 3f12.6)') 'fminsearch:   ', x_in(1), x_in(2), -fret

    ! initialize starting value for fzero
    x_in = 0.5d0

    call fzero(x_in, foc, check)

    ! print result to the screen
    write(*,'(a,3f12.6)') 'fzero:        ', x_in(1), x_in(2), -profit(x_in)

end program
