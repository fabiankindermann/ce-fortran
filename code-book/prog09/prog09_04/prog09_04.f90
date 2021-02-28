!##############################################################################
! PROGRAM StochasticGrowth2
!
! ## The stochastic growth model with an improved interpolation procedure
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
include "prog09_04m.f90"

program StochasticGrowth2

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! calculate the shock process
    call discretize_AR(rho, 0d0, sigma_eps, eta, pi)

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    do is = 1, NS
        c(:, is) = exp(eta(is))*k(:)**alpha/2
        c_new(:, is) = c(:, is)
    enddo

    ! iterate until value function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call interpolate()

        ! calculate decisions for every gridpoint and shock level
        do ik = 0, NK
            do is = 1, NS

                ! set starting value and communicate resource level
                x_in = c(ik, is)
                k_com = k(ik)
                is_com = is

                ! find the optimal consumption level
                call fzero(x_in, foc, check)
                if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                ! get optimal consumption and value function
                c_new(ik, is) = x_in

            enddo
        enddo

        ! get convergence level
        con_lev = maxval(abs(c_new(:, :) - c(:, :))/max(abs(c(:, :)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        c = c_new
    enddo

    write(*,*)'No Convergence'

contains


    ! For interpolating the RHS of the first order condition
    subroutine interpolate()

        use toolbox

        implicit none
        integer :: is, ik, is_p

        do is = 1, NS
            do ik = 0, NK

                ! calculate the RHS of the first order condition
                RHS(ik, is) = 0d0
                do is_p = 1, NS
                    RHS(ik, is) = RHS(ik, is) + pi(is, is_p)* &
                        (1d0+exp(eta(is_p))*alpha*k(ik)**(alpha-1d0)-delta)*c_new(ik, is_p)**(-1d0/gamma)
                enddo
                RHS(ik, is) = (beta*RHS(ik, is))**(-gamma)
            enddo

            ! interpolate
            call spline_interp(RHS(:, is), coeff_r(:, is))
        enddo

    end subroutine


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot)
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, err_temp

        ! end timer
        call toc()

        do is = 1, NS
            call spline_interp(c_new(:, is), coeff_c(:, is))
        enddo

        call simulate_AR(pi, is_t)

        ! calculate the time path of consumption and capital numerically
        eta_t(0) = eta(is_t(0))
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c(:, is_t(0)), k_l, k_u)
        do it = 1, TT
            eta_t(it) = eta(is_t(it))
            k_t(it) = (1d0-delta)*k_t(it-1) + exp(eta_t(it-1))*k_t(it-1)**alpha - c_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c(:, is_t(it)), k_l, k_u)
        enddo

        ! calculate plot data for value and policy function
        do is = 1, NS
            do ik = 0, nplot
                kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
                cplot(ik) = spline_eval(kplot(ik), coeff_c(:, is), k_l, k_u)
            enddo
            call plot(kplot, cplot)
        enddo
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function c(k_t)')

        ! plot time series data
        call plot((/(dble(it),it=0,TT)/), exp(eta_t))
        call execplot(xlabel='Time t', ylabel='Technology level exp(eta_t)')

        call plot((/(dble(it),it=0,TT)/), c_t)
        call execplot(xlabel='Time t', ylabel='Consumption c_t')

        call plot((/(dble(it),it=0,TT)/), k_t)
        call execplot(xlabel='Time t', ylabel='Capital k_t')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            do is = 1, NS
                is_com = is
                c_err = spline_eval(k_com, coeff_c(:, is), k_l, k_u)
                err_temp = abs(foc(c_err)/c_err)
                if(err_temp > err)err = err_temp
            enddo
        enddo
        write(*,'(a, es15.7)')'Euler equation error:',err

        ! quit program
        stop

    end subroutine

end program
