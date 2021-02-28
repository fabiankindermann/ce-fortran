!##############################################################################
! PROGRAM Ramsey
!
! ## The standard Ramsey growth model
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
include "prog09_01m.f90"

program Ramsey

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    c = k**alpha - delta*k

    ! iterate until policy function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call spline_interp(c, coeff_c)

        ! calculate decisions for every gridpoint
        do ik = 0, NK

            ! set starting value and communicate resource level
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
        endif

        c = c_new
    enddo

    write(*,*)'No Convergence'

contains


    ! For creating output plots
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
        c_stat = y_stat - delta*k_stat

        ! calculate the time path of consumption and capital numerically
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c, k_l, k_u)
        y_t(0) = k_t(0)**alpha
        do it = 1, TT
            k_t(it) = (1d0-delta)*k_t(it-1) + k_t(it-1)**alpha - c_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c, k_l, k_u)
            y_t(it) = k_t(it)**alpha
        enddo

        call plot((/(dble(it),it=0,TT)/), k_t)
        call plot((/0d0, dble(TT)/), (/k_stat, k_stat/))
        call execplot(xlabel='Time t', ylabel='Capital Stock k_t')

        call plot((/(dble(it),it=0,TT)/), c_t, legend='Consumption')
        call plot((/(dble(it),it=0,TT)/), c_t+delta*k_t, legend='Cons. + repl.Inv')
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

        ! quit program
        stop

    end subroutine

end program
