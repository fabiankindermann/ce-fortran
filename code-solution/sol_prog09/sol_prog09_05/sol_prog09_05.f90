!##############################################################################
! PROGRAM RamseyTax
!
! ## Ramsey model with variable labor supply and taxes
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog09_05m.f90"

program RamseyTax

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    l = 0.5d0
    c = max(nu/(1d0-nu)*(1d0-alpha)*(k/l)**alpha*(1d0-l), 1d-10)
    coeff_r = 0d0
    tax = 1
    tau_r_init = 0d0
    tau_w_init = 0d0
    if(tax == 1)then
        tau_r = tau_r_init
    else
        tau_w = tau_w_init
    endif

    ! compute steady state
    call steadystate()

    ! iterate until value function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call interpolate()

        ! calculate optimal decision for every gridpoint and shock level
        do ik = 0, NK

            ! set communication variables
            ik_com = ik
            k_com = k(ik)

            ! set starting value
            x_in = l(ik)

            ! find the optimal consumption level
            call fzero(x_in, foc, check)
            if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

            ! get optimal consumption and value function
            l(ik) = max(min(x_in, 1d0-1d-10), 1d-10)
            if(tax == 1)then
                tau_w(ik) = tau_w_com
            else
                tau_r(ik) = tau_r_com
            endif

            c_new(ik) = nu/(1d0-nu)*(1d0-tau_w(ik))*(1d0-alpha)*(k(ik)/l(ik))**alpha*(1d0-l(ik))

        enddo

        ! get convergence level
        con_lev = maxval(abs(c_new(:) - c(:))/max(abs(c(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        c = c_new
    enddo

    write(*,*)'No Convergence'

contains


    ! compute steady state
    subroutine steadystate()

        implicit none
        real*8 :: ck, lk, l_stat, k_stat, y_stat, c_stat, x_in(3)
        logical :: check

        ! compute steady state without taxes
        lk = ((1d0/beta-1d0+delta)/alpha)**(1d0/(1d0-alpha))
        ck = lk**(1d0-alpha)-delta
        k_stat = (ck*(1d0-nu)*lk**alpha/(nu*(1d0-alpha)) + lk)**(-1d0)
        l_stat = lk*k_stat
        c_stat = ck*k_stat
        y_stat = k_stat**alpha*l_stat**(1d0-alpha)
        gbar = gy*y_stat

        ! set initial guess
        x_in(1) = k_stat
        x_in(2) = l_stat
        x_in(3) = 0d0

        ! find root
        call fzero(x_in, foc_ss, check)
        if(check)stop 'Error: fzero did not converge'

        ! output
        c_stat = x_in(1)**alpha*x_in(2)**(1d0-alpha) - gbar - delta*x_in(1)
        V_stat = (c_stat**nu*(1d0-x_in(2))**(1d0-nu))**egam/egam/(1d0-beta)

        write(*,'(a)')'      taur      tauw         k         l         c         V'
        if(tax == 1)write(*,'(6f10.2)')tau_r_init, x_in(3), x_in(1), x_in(2), c_stat, V_stat
        if(tax == 0)write(*,'(6f10.2)')x_in(3), tau_w_init, x_in(1), x_in(2), c_stat, V_stat

        k0 = k_stat

    end subroutine


    ! for interpolating the RHS of the first order condition
    subroutine interpolate()

        implicit none

        integer :: ik
        real*8 :: marg_util

        ! calculate the RHS of the first order condition
        do ik = 0, NK
            marg_util = nu*(c(ik)**nu*(1d0-l(ik))**(1d0-nu))**egam/c(ik)
            RHS(ik) = ((beta*(1d0+(1d0-tau_r(ik))*alpha*(k(ik)/l(ik))**(alpha-1d0)-delta))*marg_util)**(-gamma)
        enddo

        ! interpolate coefficients
        call spline_interp(RHS(:), coeff_r(:))

    end subroutine


    ! For creating output plots.
    subroutine output()

        implicit none

        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot)
        real*8 :: lplot(0:nplot), tau_r_plot(0:nplot), tau_w_plot(0:nplot)
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, l_err, err_temp

        ! end timer
        call toc()

        ! interpolate policy functions
        call spline_interp(c_new(:), coeff_c(:))
        call spline_interp(l(:), coeff_l(:))
        call spline_interp(tau_r(:), coeff_tr(:))
        call spline_interp(tau_w(:), coeff_tw(:))

        ! calculate plot data for consumption policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            cplot(ik) = spline_eval(kplot(ik), coeff_c(:), k_l, k_u)
        enddo

        call plot(kplot, cplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function c(k_t)')

        ! calculate plot data for labor policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            lplot(ik) = spline_eval(kplot(ik), coeff_l(:), k_l, k_u)
        enddo
        call plot(kplot, lplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function l(k_t)')

        do ik = 0, nplot
            tau_r_plot(ik) = spline_eval(kplot(ik), coeff_tr(:), k_l, k_u)
            tau_w_plot(ik) = spline_eval(kplot(ik), coeff_tw(:), k_l, k_u)
        enddo
        call plot(kplot, tau_r_plot, legend='tau_r(k)')
        call plot(kplot, tau_w_plot, legend='tau_w(k)')
        call execplot(xlabel='Capital', ylabel='Taxes')

        ! calculate the time path of consumption and capital numerically
        k_t(0) = k0
        tau_w_t(0) = spline_eval(k_t(0), coeff_tw, k_l, k_u)
        tau_r_t(0) = spline_eval(k_t(0), coeff_tr, k_l, k_u)
        c_t(0) = spline_eval(k_t(0), coeff_c, k_l, k_u)
        l_t(0) = spline_eval(k_t(0), coeff_l, k_l, k_u)
        y_t(0) = k_t(0)**alpha*l_t(0)**(1d0-alpha)
        V_trans = (c_t(0)**nu*(1d0-l_t(0))**(1d0-nu))**egam/egam

        do it = 1, TT
            k_t(it) = (1d0-delta)*k_t(it-1) + y_t(it-1) - c_t(it-1) - gbar
            tau_w_t(it) = spline_eval(k_t(it), coeff_tw, k_l, k_u)
            tau_r_t(it) = spline_eval(k_t(it), coeff_tr, k_l, k_u)
            c_t(it) = spline_eval(k_t(it), coeff_c, k_l, k_u)
            l_t(it) = spline_eval(k_t(it), coeff_l, k_l, k_u)
            y_t(it) = k_t(it)**alpha*l_t(it)**(1d0-alpha)
            V_trans = V_trans + beta**it*(c_t(it)**nu*(1d0-l_t(it))**(1d0-nu))**egam/egam
        enddo
        call plot((/(dble(it),it=0,TT)/), l_t, legend='l')
        call plot((/(dble(it),it=0,TT)/), c_t, legend='c')
        call execplot(xlabel='Time', ylabel='Labor, Consumption')

        call plot((/(dble(it),it=0,TT)/), k_t, legend='k')
        call plot((/(dble(it),it=0,TT)/), y_t, legend='y')
        call execplot(xlabel='Time', ylabel='k, y')

        call plot((/(dble(it),it=0,TT)/), tau_r_t, legend='tau_r')
        call plot((/(dble(it),it=0,TT)/), tau_w_t, legend='tau_w')
        call execplot(xlabel='Time', ylabel='Labor, Consumption')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            c_err = spline_eval(k_com, coeff_c(:), k_l, k_u)
            l_err = spline_eval(k_com, coeff_l(:), k_l, k_u)
            err_temp = abs(foc(l_err)/c_err)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(/a, es15.7)')'Euler equation error:',err

        ! quit program
        stop

    end subroutine

end program
