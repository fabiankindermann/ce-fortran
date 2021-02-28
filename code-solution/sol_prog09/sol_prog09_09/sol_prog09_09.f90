!##############################################################################
! PROGRAM Hansen
!
! ## Stochastic growth model with indivisible labor - Hansen(1985)
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
include "sol_prog09_09m.f90"

program Hansen

    ! modules
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
        kappa(:, is) = 0.5d0
        c(:, is) = exp(eta(is))*k(:)**alpha*(kappa(:, is)*h0)**(1d0-alpha)
        c_new(:, is) = c(:, is)
    enddo

    ! iterate until value function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call interpolate()

        ! calculate optimal decision for every gridpoint and shock level
        do ik = 0, NK
            do is = 1, NS

                ! initialize starting value and communicate resource level
                x_in = kappa(ik, is)
                k_com = k(ik)
                is_com = is

                ! find the optimal consumption level
                call fzero(x_in, foc, check)
                if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                ! get optimal consumption and value function
                kappa(ik, is) = x_in
                c_new(ik, is) = -(1d0-alpha)*exp(eta(is))*(k(ik)/(x_in*h0))**alpha*h0/(A*log(1d0-h0))

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
                                  (1d0+alpha*exp(eta(is_p))*(k(ik)/(kappa(ik, is_p)*h0))&
                                  **(alpha-1d0) - delta)/c_new(ik, is_p)
                enddo
                RHS(ik, is) = beta*RHS(ik, is)
            enddo

            ! interpolate
            call spline_interp(RHS(:, is), coeff_r(:, is))
        enddo

    end subroutine


    ! For creating output plots.
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), cplot(0:nplot), kappaplot(0:nplot)
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, kappa_err, err_temp

        ! end timer
        call toc()

        ! interpolate policy functions c(k,eta) and alpha(k,eta)
        do is = 1, NS
            call spline_interp(c_new(:, is), coeff_c(:, is))
            call spline_interp(kappa(:, is), coeff_kappa(:, is))
        enddo

        ! Analytical Steady State, the first two equations are provided in the paper
        hbar = ((1d0-alpha)*(1d0/beta-1d0+delta))/(((-A*(log(1d0-h0)))/h0)*(1d0/beta-1d0+delta-alpha*delta))
        kbar = ((1d0/beta-1d0+delta)/alpha)**(1d0/(alpha-1d0))*hbar
        ybar = exp(eta_t(TT))*kbar**alpha*hbar**(1d0-alpha)
        cbar = -((exp(eta_t(TT))*(kbar/(hbar))**alpha*(1d0-alpha)*h0)/(A*log(1d0-h0)))
        ibar = ybar-cbar
        rbar = exp(eta_t(TT))*alpha*(hbar/kbar)**(1d0-alpha)-delta
        wbar = exp(eta_t(TT))*(1d0-alpha)*(kbar/hbar)**alpha*h0

        ! compare the long run stochastic steady state with the analytical one.
        is_t = floor(dble(NS)/2d0)+1
        call simulate_economy
        k0 = k_t(TT)

        write(*,'(a)')'                              c       i       k       h       r       w'
        write(*,'(a,7f8.2)')'Long Run Steady State:  ', c_t(TT)/y_t(TT)*100d0, i_t(TT)/y_t(TT)*100d0, &
            k_t(TT)/(4d0*y_t(TT))*100d0, kappa_t(TT)*h0*100d0, ((1d0+r_t(TT))**4-1d0)*100d0, w_t(TT)

        write(*,'(a,7f8.2)')'Analytical Steady State:', cbar/ybar*100d0, ibar/ybar*100d0, &
            (kbar/(4d0*ybar))*100d0, hbar*100d0, ((1d0+rbar)**4-1d0)*100d0, wbar

        ! simulate a series of shocks and then determine normalized time series
        call simulate_AR(pi, is_t)

        ! simulate the economy forward
        call simulate_economy()

        ! calculate business cycle statistics
        call bc_statistics()

        ! write output
        write(*,'(//a/)')'BUSINESS CYCLE MOMENTS'
        write(*,'(a)')'                               c       i       k   kappa       r       w'
        write(*,'(a,6f8.2)')'Mean:                   ',(/E_c, E_i, E_k/)/E_y*100d0, (/E_kappa, E_r/)*100d0, E_w

        write(*,'(a,6f8.2)')'CV (in % of CV of GDP): ',(/CV_c, CV_i, CV_k, CV_kappa, CV_r, CV_w/)/CV_y
        write(*,'(a,6f8.2)')'Corr. with GDP:         ',Cor_c, Cor_i, Cor_k, Cor_kappa, Cor_r, Cor_w

        ! plot some normalized time series
        call plot((/(dble(it),it=0,5000)/), c_t(0:5000))
        call execplot(xlabel='Time t', ylabel='Consumption c_t')

        call plot((/(dble(it),it=0,5000)/), kappa_t(0:5000)*h0)
        call execplot(xlabel='Time t', ylabel='Labor Supply h_t')

        call plot((/(dble(it),it=0,5000)/), k_t(0:5000))
        call execplot(xlabel='Time t', ylabel='Capital k_t')

        ! calculate plot data for consumption policy function
        do is = 1, NS
            do ik = 0, nplot
                kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
                cplot(ik) = spline_eval(kplot(ik), coeff_c(:, is), k_l, k_u)
            enddo
            call plot(kplot, cplot)
        enddo
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function c(k_t)')

        ! calculate plot data for labor policy function
        do is = 1, NS
            do ik = 0, nplot
                kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
                kappaplot(ik) = spline_eval(kplot(ik), coeff_kappa(:, is), k_l, k_u)
            enddo
            call plot(kplot, kappaplot)
        enddo
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function kappa(k_t)')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            do is = 1, NS
                is_com = is
                c_err = spline_eval(k_com, coeff_c(:, is), k_l, k_u)
                kappa_err = spline_eval(k_com, coeff_kappa(:, is), k_l, k_u)
                err_temp = abs(foc(kappa_err)/c_err)
                if(err_temp > err)err = err_temp
            enddo
        enddo
        write(*,'(/a, es15.7)')'Euler equation error:', err

        ! quit program
        stop

    end subroutine


    ! simulates the economy for a given series of shocks
    subroutine simulate_economy()

        implicit none

        ! calculate the time path of all economic variables
        eta_t(0) = eta(is_t(0))
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c(:, is_t(0)), k_l, k_u)
        kappa_t(0) = spline_eval(k_t(0), coeff_kappa(:, is_t(0)), k_l, k_u)
        y_t(0) = exp(eta_t(0))*k_t(0)**alpha*(kappa_t(0)*h0)**(1d0-alpha)
        i_t(0) = y_t(0) - c_t(0)
        r_t(0) = exp(eta_t(0))*alpha*((kappa_t(0)*h0)/k_t(0))**(1d0-alpha) - delta
        w_t(0) = exp(eta_t(0))*(1d0-alpha)*(k_t(0)/(kappa_t(0)*h0))**alpha*h0
        do it = 1, TT
            eta_t(it) = eta(is_t(it))
            k_t(it) = (1d0-delta)*k_t(it-1) + i_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c(:, is_t(it)), k_l, k_u)
            kappa_t(it) = spline_eval(k_t(it), coeff_kappa(:, is_t(it)), k_l, k_u)
            y_t(it) = exp(eta_t(it))*k_t(it)**alpha*(kappa_t(it)*h0)**(1d0-alpha)
            i_t(it) = y_t(it)-c_t(it)
            r_t(it) = exp(eta_t(it))*alpha*((kappa_t(it)*h0)/k_t(it))**(1d0-alpha) - delta
            w_t(it) = exp(eta_t(it))*(1d0 - alpha)*(k_t(it)/(kappa_t(it)*h0))**alpha*h0
        enddo

    end subroutine


    ! calculates business cycle statistics
    subroutine bc_statistics()

        implicit none

        ! calculate expectations
        E_c = sum(c_t)/dble(TT+1)
        E_i = sum(i_t)/dble(TT+1)
        E_k = sum(k_t)/dble(TT+1)
        E_kappa = sum(kappa_t)/dble(TT+1)
        E_r = sum(r_t)/dble(TT+1)
        E_w = sum(w_t)/dble(TT+1)
        E_y = sum(y_t)/dble(TT+1)

        ! normalization on annual level
        E_k = E_k/4d0
        E_r = (1d0+E_r)**4-1d0

        ! calculate coefficient of variation
        CV_c = CofVar(c_t)
        CV_i = CofVar(i_t)
        CV_k = CofVar(k_t)
        CV_kappa = CofVar(kappa_t)
        CV_r = CofVar(r_t)
        CV_w = CofVar(w_t)
        CV_y = CofVar(y_t)

        ! calculate correlation coefficients of logs
        Cor_c = correlation(c_t, y_t)
        Cor_i = correlation(i_t, y_t)
        Cor_k = correlation(k_t, y_t)
        Cor_kappa = correlation(kappa_t, y_t)
        Cor_r = correlation(r_t, y_t)
        Cor_w = correlation(w_t, y_t)

    end subroutine


    ! Calculates coefficient of variation
    function CofVar(x)

        implicit none
        real*8, intent(in) :: x(0:TT)
        real*8 :: CofVar, mean
        integer :: it

        mean = sum(x)/dble(TT+1)
        CofVar = 0d0
        do it = 0, TT
            CofVar = CofVar + (x(it)-mean)**2
        enddo
        CofVar = sqrt(CofVar/dble(TT))/max(mean, 1d-10)

    end function


    ! Calculates correlation coefficient
    function correlation(x, y)

        implicit none
        real*8, intent(in) :: x(0:TT), y(0:TT)
        real*8 :: correlation, variance1, variance2, mean1, mean2
        integer :: it

        mean1 = sum(x)/dble(TT+1)
        mean2 = sum(y)/dble(TT+1)
        variance1   = 0d0
        variance2   = 0d0
        correlation = 0d0
        do it = 0, TT
            variance1   = variance1   + (x(it)-mean1)**2
            variance2   = variance2   + (y(it)-mean2)**2
            correlation = correlation + (x(it)-mean1)*(y(it)-mean2)
        enddo
        variance1   = variance1/dble(TT)
        variance2   = variance2/dble(TT)
        correlation = correlation/dble(TT)/max(sqrt(variance1*variance2), 1d-10)

    end function

end program
