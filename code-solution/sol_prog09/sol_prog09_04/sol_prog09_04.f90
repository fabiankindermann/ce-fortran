!##############################################################################
! PROGRAM RamseyLabor
!
! ## Ramsey model with variable labor supply and value function iteration
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog09_04m.f90"

program RamseyLabor

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize grid, policy and value function
    call grid_Cons_Equi(k, k_l, k_u)
    l = 0.5d0
    c = k**alpha*l**(1d0 - alpha) - delta*k
    V = 0d0
    coeff_V = 0d0

    ! compute steady state
    call steadystate()

    ! iterate until value function converges
    do iter = 1, itermax

        ! calculate optimal decision for every gridpoint
        do ik = 0, NK

            ! initialize starting value (labor) and communicate resource level
            x_in = l(ik)
            k_com = k(ik)

            ! find the optimal labour level
            call fminsearch(x_in, fret, 1d-10, 1d0-1d-10, utility)

            ! get optimal consumption and value function
            c(ik) = (nu/(1d0-nu))*(1d0-x_in)*(1d0-alpha)*(k_com/x_in)**alpha
            V_new(ik) = -fret
            l(ik) = x_in

        enddo

        ! interpolate coefficients
        call spline_interp((egam*V_new)**(1d0/egam), coeff_V)

        ! get convergence level
        con_lev = maxval(abs(V_new(:) - V(:))/max(abs(V(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        V = V_new

    enddo

    write(*,*)'No Convergence'

contains


    ! compute steady state
    subroutine steadystate()

        implicit none
        real*8 :: x_in(2), ck, lk
        logical :: check

        x_in(1) = 4d0
        x_in(2) = 0.1d0

        ! find root
        call fzero(x_in, focs, check)
        if(check)stop 'Error: fzero did not converge'

        ! output
        write(*,'(a,f10.4)')'Capital: ', x_in(1)
        write(*,'(a,f10.4)')'Labor:   ', x_in(2)
        write(*,'(a,f10.4)')'Consumption: ', x_in(1)**alpha*x_in(2)**(1d0-alpha) - delta*x_in(1)
        write(*,'(a,f10.4)')'Output: ', x_in(1)**alpha*x_in(2)**(1d0-alpha)

        ! analytical solution
        lk = ((1d0/beta-1d0+delta)/alpha)**(1d0/(1d0-alpha))
        ck = lk**(1d0-alpha) - delta
        k_stat = (ck*(1d0-nu)*lk**alpha/(nu*(1d0-alpha)) + lk)**(-1d0)
        l_stat = lk*k_stat
        c_stat = ck*k_stat
        y_stat = k_stat**alpha*l_stat**(1d0-alpha)

        write(*,'(/a,f10.4)')'Capital: ', k_stat
        write(*,'(a,f10.4)')'Labor:   ', l_stat
        write(*,'(a,f10.4)')'Consumption: ', c_stat
        write(*,'(a,f10.4)')'Output: ', y_stat

    end subroutine


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot), lplot(0:nplot)

        ! end timer
        call toc()

        ! interpolate policy function
        call spline_interp(c, coeff_c)
        call spline_interp(l, coeff_l)

        ! calculate the time path of consumption and capital numerically
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c, k_l, k_u)
        l_t(0) = spline_eval(k_t(0), coeff_l, k_l, k_u)
        y_t(0) = k_t(0)**alpha
        do it = 1, TT
            k_t(it) = (1d0-delta)*k_t(it-1) + k_t(it-1)**alpha*l_t(it-1)**(1d0-alpha) - c_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c, k_l, k_u)
            l_t(it) = spline_eval(k_t(it), coeff_l, k_l, k_u)
            y_t(it) = k_t(it)**alpha*l_t(it)**(1d0-alpha)
        enddo

        call plot((/(dble(it),it=0,TT)/), k_t)
        call plot((/0d0, dble(TT)/), (/k_stat, k_stat/))
        call execplot(xlabel='Time t', ylabel='Capital Stock k_t')

        call plot((/(dble(it),it=0,TT)/), c_t, legend='Consumption')
        call plot((/(dble(it),it=0,TT)/), l_t, legend='Labour')
        call plot((/(dble(it),it=0,TT)/), c_t+delta*k_t, legend='Cons. + repl.Inv')
        call plot((/(dble(it),it=0,TT)/), y_t, legend='Output')
        call execplot(xlabel='Time t', ylabel='Goods Market')

        ! calculate plot data for consumption policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            cplot(ik) = spline_eval(kplot(ik), coeff_c(:), k_l, k_u)
        enddo

        ! plot policy function
        call plot(kplot, cplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function c(k_t)')

        ! calculate plot data for labor policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            lplot(ik) = spline_eval(kplot(ik), coeff_l, k_l, k_u)
        enddo

        ! plot policy function
        call plot(kplot, lplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function l(k_t)')

        ! quit program
        stop

    end subroutine

end program
