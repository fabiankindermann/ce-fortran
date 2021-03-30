!##############################################################################
! PROGRAM MarketSplines
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
include "sol_prog02_09m.f90"

program MarketSplines

    use globals
    use toolbox

    implicit none

    ! declaration of variables
    integer :: ip
    real*8 :: q_in
    logical :: check

    ! initialize coefficients for interpolation to zero
    coeff_q = 0d0

    ! initialize grid for interpolation
    call grid_Cons_Equi(p, p_l, p_u)

    ! calculate quantities q at each gridpoint
    do ip = 0, N

        ! set initial guess
        q_in = 1d0

        ! set communicatation variable
        ip_com = ip

        ! call fzero
        call fzero(q_in, foc, check)
        if(check)stop 'Error: fzero did not converge'

        ! copy solution at grid point
        q(ip) = q_in

    enddo

    ! print values for each gridpoint q(ip)
    write(*,'(a)')'Individual quantities at each gridpoint:'
    write(*,'(a)')' Price  |    Quantity'
    do ip = 0, N
        write(*,'(f7.4,a,f13.10)')p(ip), ' | ', q(ip)
    enddo

    ! compute coefficients for spline
    call spline_interp(q, coeff_q)

    ! set up grid for plot
    call grid_Cons_Equi(p_plot, p_l, p_u)

    ! interpolate supply for plot
    q_s_plot = (/(spline_eval(p_plot(ip), coeff_q, p_l, p_u), ip = 0, NP)/)

    ! compute demand
    q_d_plot = p_plot**(-eta)

    ! initialize plot
    call plot(p_plot, q_d_plot, legend='Market Demand')
    call plot(p_plot, q_s_plot, legend='Individual Supply')
    call plot(p_plot, m*q_s_plot, legend='Aggregate Supply')

    ! execute plot
    call execplot(xlim=(/0.1d0, 3d0/), ylim=(/0d0, 2d0/), &
                  xlabel='Price', ylabel='Demand/Quantity')

    ! initial guess to compute the equilibrium price
    p_equi = 1d0

    ! call fzero
    call fzero(p_equi, market, check)
    if(check)stop 'Error: fzero did not converge'

    ! print output
    write(*,'(a,f13.10)')'Equilibrium Price:', p_equi
    write(*,'(a,f13.10)')'Equilibrium Quantity:', p_equi**(-eta)

end program
