!##############################################################################
! PROGRAM StochasticGrowthVFIEP
!
! ## Stochastic growth model with value function iteration and end. gridpoints
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog09_07m.f90"

program StochasticGrowthVFIEP

    use globals
    use toolbox

    implicit none

    ! start timer for value function iteration
    call tic()

    ! calculate the shock process
    call discretize_AR(rho, 0d0, sigma_eps, eta, pi)

    ! initialize grid on capital tomorrow
    call grid_Cons_Equi(k, k_l, k_u)
    do is = 1, NS
        c(:, is) = 0.5d0*exp(eta(is))*k(:)**alpha/2
    enddo
    V = 1d0
    coeff_EV = 0.1d0

    ! iterate until value function converges
    do iter = 1, itermax

        ! calculate optimal decision for every gridpoint
        do ik = 0, NK
            do is = 1, NS

                ! initialize starting value and communicate cash-on-hand
                X_com = (1d0-delta)*k(ik) + exp(eta(is))*k(ik)**alpha
                x_in = X_com - c(ik, is)
                is_com = is

                ! find the optimal future capital
                call fminsearch(x_in, fret, k_l, min(X_com, k_u), utility)

                ! get optimal consumption function
                c(ik, is) = X_com - x_in
                V_new(ik, is) = -fret

            enddo
        enddo

        ! interpolate coefficients
        call interpolate()

        ! get convergence level
        con_lev = maxval(abs(V_new(:, :) - V(:, :))/max(abs(V(:, :)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
            exit
        endif

        V = V_new
    enddo

    if(con_lev > sig)write(*,*)'No Convergence with value function iteration'

    ! initialize cash-on-hand for tomorrow
    do is = 1, NS
        X(:, is) = exp(eta(is))*k(:)**alpha + (1d0-delta)*k(:)
    enddo

    ! start timer for value function iteration
    call tic()

    ! initialize value functions
    do is = 1, NS
        VX(:, is) = (X(:, is) - k(:))**egam/egam
    enddo

    V = 0d0
    do is = 1, NS
        do is_p = 1, NS
            V(:, is) = V(:, is) + pi(is, is_p)*VX(:, is_p)
        enddo
    enddo

    ! iterate until value function converges
    do iter = 1, itermax

        ! compute derivative of value function
        do ik = 0, NK
            do is = 1, NS
                if(ik == 0)then
                    dV(ik, is) = (V(ik+1, is) - V(ik, is))/(k(ik+1) - k(ik))
                elseif(ik == NK)then
                    dV(ik, is) = (V(ik, is) - V(ik-1, is))/(k(ik) - k(ik-1))
                else
                    dV(ik, is) = 0.5d0*((V(ik+1, is) - V(ik, is))/(k(ik+1) - k(ik)) + &
                                        (V(ik, is) - V(ik-1, is))/(k(ik) - k(ik-1)))
                endif
            enddo
        enddo

        ! determine endogenous consumption function, cash-on-hand and respective value function
        c_endog = 0d0
        do ik = 0, NK
            do is = 1, NS
                do is_p = 1, NS
                    c_endog(ik, is) = c_endog(ik, is) + pi(is, is_p)*dV(ik, is_p)
                enddo
                c_endog(ik, is) = (beta*c_endog(ik, is))**(-gamma)
                X_endog(ik, is) = c_endog(ik, is) + k(ik)
                VX_endog(ik, is) = c_endog(ik, is)**egam/egam + beta*V(ik, is)
            enddo
        enddo

        ! update value function on resources by interpolation
        do ik = 0, NK
            do is = 1, NS
                VX(ik, is) = linint_Gen(X(ik, is), X_endog(:, is), VX_endog(:, is), ik)
            enddo
        enddo

        ! update initial value function on capital
        V_new = 0d0
        do is = 1, NS
            do is_p = 1, NS
                V_new(:, is) = V_new(:, is) + pi(is, is_p)*VX(:, is_p)
            enddo
        enddo

        ! get convergence level
        con_lev = maxval(abs(V_new(:, :) - V(:, :))/max(abs(V(:, :)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then

            ! calculate endogenous capital for every shock level
            do ik = 0, NK
                do is = 1, NS

                    ! set starting value and communicate resource level
                    x_in = k(ik)
                    X_com = X_endog(ik, is)
                    is_com = is

                    call fzero(x_in, fkend, check)
                    if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                    ! get endogenous capital stock
                    k_endog(ik, is) = x_in

                enddo
            enddo

            ! compute policy function by interpolation
            do ik = 0, NK
                do is = 1, NS
                    c(ik, is) = linint_Gen(k(ik), k_endog(:, is), c_endog(:, is), ik)
                enddo
            enddo

            call output()
            stop
        endif

        V = V_new
    enddo

    if(con_lev > sig)write(*,*)'No Convergence with endogenous gridpoints method'

contains


    ! For interpolating the expected value function
    subroutine interpolate()

        implicit none
        integer :: is, ik, is_p

        do is = 1, NS
            do ik = 0, NK

                ! calculate the expected value function
                EV(ik, is) = 0d0
                do is_p = 1, NS
                    EV(ik, is) = EV(ik, is) + pi(is, is_p)*V_new(ik, is_p)
                enddo
            enddo

            ! interpolate
            call spline_interp((egam*EV(:, is))**(1d0/egam), coeff_EV(:, is))
        enddo

    end subroutine


    ! For creating output plots.
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot)
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, err_temp

        ! end timer
        call toc()

        ! simulate a series of shocks
        call simulate_AR(pi, is_t)

        ! interpolate coefficients
        do is = 1, NS
            call spline_interp(c(:, is), coeff_c(:, is))
        enddo

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
                X_com = (1d0-delta)*k_com + exp(eta(is))*k_com**alpha
                is_com = is
                c_err = spline_eval(k_com, coeff_c(:, is), k_l, k_u)
                err_temp = abs(foc(c_err)/c_err)
                if(err_temp > err)err = err_temp
            enddo
        enddo
        write(*,'(a, es15.7)')'Euler equation error:',err

    end subroutine

end program
