!##############################################################################
! PROGRAM RBC_cost
!
! ## Calculates the cost of business cycle fluctuations
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
include "prog09_06m.f90"

program RBC_cost

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
    do is = 1, ns
        l(:, is) = 0.3d0
        c(:, is) = exp(eta(is))*k(:)**alpha*l(:, is)**(1d0-alpha)
        c_new(:, is) = c(:, is)
    enddo

    ! iterate until policy function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call interpolate()

        ! calculate decisions for every gridpoint and shock level
        do ik = 0, NK
            do is = 1, NS

                ! sets starting value and communicate resource level
                x_in = l(ik, is)
                k_com = k(ik)
                is_com = is

                ! find the optimal labor supply level
                call fzero(x_in, foc, check)
                if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                ! get optimal consumption function
                l(ik, is) = x_in
                c_new(ik, is) = nu/(1d0-nu)*(1d0-x_in)*exp(eta(is))*(1d0-alpha)*(k(ik)/x_in)**alpha

            enddo
        enddo

        ! get convergence level
        con_lev = maxval(abs(c_new(:, :) - c(:, :))/max(abs(c(:, :)), 1d-10))
        if(mod(iter,50) == 0)write(*,'(i5,2x,f20.7)')iter, con_lev

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
        real*8 :: marg_util

        do is = 1, NS
            do ik = 0, NK

                ! calculate the RHS of the first order condition
                RHS(ik, is) = 0d0
                do is_p = 1, NS
                    marg_util = nu*(c_new(ik, is_p)**nu*(1d0-l(ik, is_p))**(1d0-nu))**egam/c_new(ik, is_p)
                    RHS(ik, is) = RHS(ik, is) + pi(is, is_p)* &
                        (1d0+exp(eta(is_p))*alpha*(k(ik)/l(ik, is_p))**(alpha-1d0)-delta)*marg_util
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
        real*8 :: Omega, k_star, c_star, l_star, u_star
        real*8 :: CV_leis, cor_leis
        real*8 :: u_mu, u_hat, bar_nu
        real*8 :: bc_lev, bc_c, bc_leis, bc_cor, bc_tot
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, l_err, err_temp

        ! end timer
        call toc()

        ! interpolate policy functions
        do is = 1, NS
            call spline_interp(c_new(:, is), coeff_c(:, is))
            call spline_interp(l(:, is), coeff_l(:, is))
        enddo

        ! determine the stochastic steady state
        is_t = floor(dble(NS)/2d0)+1
        call simulate_economy()
        k0 = k_t(TT)

        ! simulate a series of shocks
        call simulate_AR(pi, is_t)

        ! simulate the economy forward
        call simulate_economy()

        ! calculate business cycle statistics
        call bc_statistics()

        ! write output
        write(*,'(//a/)')'BUSINESS CYCLE MOMENTS'
        write(*,'(a)')'                               c       i       k       l       r       w'
        write(*,'(a,6f8.2)')'Mean:                   ',(/E_c, E_i, E_k/)/E_y*100d0, (/E_l, E_r/)*100d0, E_w

        write(*,'(a,6f8.2)')'CV (in % of CV of GDP): ',(/CV_c, CV_i, CV_k, CV_l, CV_r, CV_w/)/CV_y
        write(*,'(a,6f8.2)')'Corr. with GDP:         ',Cor_c, Cor_i, Cor_k, Cor_l, Cor_r, Cor_w


        !########################################################
        ! GET THE COST OF THE BUSINESS CYCLE
        !########################################################

        ! get CV and correlation
        CV_leis  = CofVar(1d0-l_t)
        Cor_leis = correlation(c_t, 1d0-l_t)

        ! determine steady state values under certainty
        Omega  = exp(0.5d0*sigma_eps/(1d0-rho**2))
        k_star = ((Omega*alpha)/(1d0/beta-1d0+delta))**(1d0/(1d0-alpha))
        l_star = nu*Omega*(1d0-alpha)*k_star**alpha/ &
            ((1d0-nu*alpha)*Omega*k_star**alpha - (1d0-nu)*delta*k_star)
        c_star = (Omega*k_star**alpha-delta*k_star)*l_star

        ! get different utility levels
        u_star = (c_star**nu*(1d0-l_star)**(1d0-nu))**egam/egam
        u_mu   = (E_c**nu*(1d0-E_l)**(1d0-nu))**egam/egam
        u_hat  = u_mu/u_star
        E_V    = sum((c_t**nu*(1d0-l_t)**(1d0-nu))**egam/egam)/dble(TT+1)

        ! get different components of business cycle costs
        bar_nu  = (1d0-nu)*nu*gamma
        bc_lev  = (u_hat-1d0)/(nu*egam)
        bc_c    = -u_hat/(2d0*nu*gamma)*(nu**2 + bar_nu)*CV_c**2
        bc_leis = -u_hat/(2d0*nu*gamma)*((1d0-nu)**2 + bar_nu)*CV_leis**2
        bc_cor  = -u_hat/(nu*gamma)*(nu*(1d0-nu) - bar_nu)*Cor_leis*CV_c*CV_leis
        bc_tot  = bc_lev + bc_c + bc_leis + bc_cor

        write(*,'(//a/)')'COST OF THE BUSINESS CYCLE'
        write(*,'(a)')'                     LEVEL       CV(c)       CV(l)    RHO(c,l)       TOTAL'
        write(*,'(a,12x,3f12.5)')'Values:       ',CV_c, CV_leis, Cor_leis
        write(*,'(a,5f12.5)')'Aprox. cost:  ',(/bc_lev, bc_c, bc_leis, bc_cor, bc_tot/)*100d0
        write(*,'(a,48x,f12.5//)')'Model sim.:   ',((E_V/u_star)**(1d0/(nu*egam))-1d0)*100d0

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            do is = 1, NS
                is_com = is
                c_err = spline_eval(k_com, coeff_c(:, is), k_l, k_u)
                l_err = spline_eval(k_com, coeff_l(:, is), k_l, k_u)
                err_temp = abs(foc2(l_err)/c_err)
                if(err_temp > err)err = err_temp
            enddo
        enddo
        write(*,'(a, es15.7)')'Euler equation error:',err

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
        l_t(0) = spline_eval(k_t(0), coeff_l(:, is_t(0)), k_l, k_u)
        y_t(0) = exp(eta_t(0))*k_t(0)**alpha*l_t(0)**(1d0-alpha)
        i_t(0) = y_t(0)-c_t(0)
        r_t(0) = exp(eta_t(0))*alpha*(l_t(0)/k_t(0))**(1d0-alpha) - delta
        w_t(0) = exp(eta_t(0))*(1d0-alpha)*(k_t(0)/l_t(0))**(alpha)
        do it = 1, TT
            eta_t(it) = eta(is_t(it))
            k_t(it) = (1d0-delta)*k_t(it-1) + i_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c(:, is_t(it)), k_l, k_u)
            l_t(it) = spline_eval(k_t(it), coeff_l(:, is_t(it)), k_l, k_u)
            y_t(it) = exp(eta_t(it))*k_t(it)**alpha*l_t(it)**(1d0-alpha)
            i_t(it) = y_t(it)-c_t(it)
            r_t(it) = exp(eta_t(it))*alpha*(l_t(it)/k_t(it))**(1d0-alpha) - delta
            w_t(it) = exp(eta_t(it))*(1d0-alpha)*(k_t(it)/l_t(it))**alpha
        enddo

    end subroutine


    ! calculates business cycle statistics
    subroutine bc_statistics()

        implicit none

        ! calculate expectations
        E_c = sum(c_t)/dble(TT+1)
        E_i = sum(i_t)/dble(TT+1)
        E_k = sum(k_t)/dble(TT+1)
        E_l = sum(l_t)/dble(TT+1)
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
        CV_l = CofVar(l_t)
        CV_r = CofVar(r_t)
        CV_w = CofVar(w_t)
        CV_y = CofVar(y_t)

        ! calculate correlation coefficients of logs
        Cor_c = correlation(c_t, y_t)
        Cor_i = correlation(i_t, y_t)
        Cor_k = correlation(k_t, y_t)
        Cor_l = correlation(l_t, y_t)
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
