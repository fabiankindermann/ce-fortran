!##############################################################################
! PROGRAM RamseyPop
!
! ## Ramsey model with population growth, policy and value function iteration
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog09_01m.f90"

program RamseyPop

    use globals
    use toolbox

    implicit none

    ! start timer for policy function iteration
    call tic()

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    c = k**alpha - (delta+n_p)*k

    ! iterate until policy function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call spline_interp(c, coeff_c)

        ! calculate optimal decision for every gridpoint
        do ik = 0, NK

            ! initialize starting value and communicate resource level
            x_in = c(ik)
            k_com = k(ik)

            ! find the optimal consumption level
            call fzero(x_in, foc, check)
            if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

            ! get optimal consumption and value function
            c_new(ik) = x_in

        enddo

        ! get convergence level
        con_lev = maxval(abs(c_new(:) - c(:))/max(abs(c(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
            exit
        endif

        c = c_new
    enddo

    if(con_lev > sig)write(*,*)'No Convergence with policy function iteration'

    ! start timer for value function iteration
    call tic()

    ! initialize spline, policy and value function
    c = k**alpha - (delta+n_p)*k
    V = 0d0
    coeff_v = 0d0

    ! iterate until value function converges
    do iter = 1, itermax

        ! calculate optimal decision for every gridpoint
        do ik = 0, NK

            ! initialize starting value and communicate resource level
            x_in = ((1d0-delta)*k(ik) + k(ik)**alpha - c(ik))/(1d0+n_p)
            k_com = k(ik)

            call fminsearch(x_in, fret, k_l, min((1d0 - delta)*k(ik) + &
                                                 k(ik)**alpha, k_u), utility)

            ! get optimal consumption and value function
            c(ik) = (1d0-delta)*k(ik) + k(ik)**alpha - (1d0+n_p)*x_in
            V_new(ik) = -fret

        enddo

        ! interpolate coefficients
        call spline_interp((egam*V_new)**(1d0/egam), coeff_V)

        ! get convergence level
        con_lev = maxval(abs(V_new(:) - V(:))/max(abs(V(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
            stop
        endif

        V = V_new
    enddo

    if(con_lev > sig)write(*,*)'No Convergence with value function iteration'

contains


    ! For creating output plots.
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot)
        real*8 :: theta, c_stat, k_stat, y_stat
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, err_temp

        ! end timer
        call toc()

        ! calculate stationary equilibrium values
        theta = 1d0/beta-1d0
        k_stat = (alpha/(delta+theta))**(1d0/(1d0-alpha))
        y_stat = k_stat**alpha
        c_stat = y_stat - (delta+n_p)*k_stat

        ! calculate the time path of consumption and capital numerically
        call spline_interp(c, coeff_c)

        ! calculate the time path of consumption and capital numerically
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c, k_l, k_u)
        y_t(0) = k_t(0)**alpha
        do it = 1, TT
            k_t(it) = ((1d0-delta)*k_t(it-1) + y_t(it-1) - c_t(it-1))/(1d0 + n_p)
            c_t(it) = spline_eval(k_t(it), coeff_c, k_l, k_u)
            y_t(it) = k_t(it)**alpha
        enddo

        call plot((/(dble(it),it=0,TT)/), k_t)
        call plot((/0d0, dble(TT)/), (/k_stat, k_stat/))
        call execplot(xlabel='Time t', ylabel='Capital Stock k_t')

        call plot((/(dble(it),it=0,TT)/), c_t, legend='Consumption')
        call plot((/(dble(it),it=0,TT)/), c_t+(delta+n_p)*k_t, legend='Cons. + repl.Inv')
        call plot((/(dble(it),it=0,TT)/), y_t, legend='Output')
        call execplot(xlabel='Time t', ylabel='Goods Market')

        ! calculate plot data for value and policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            cplot(ik) = spline_eval(kplot(ik), coeff_c, k_l, k_u)
        enddo

        ! plot policy function
        call plot(kplot, cplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function c(k_t)')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            c_err = spline_eval(k_com, coeff_c, k_l, k_u)
            err_temp = abs(foc(c_err)/c_err)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(a, es15.7)')'Euler equation error:',err

    end subroutine

end program
