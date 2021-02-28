!##############################################################################
! PROGRAM RBC_government
!
! ## Pro-cyclical versus constant government expenditure in the RBC model
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
include "prog09_07m.f90"

program RBC_government

    use globals
    use toolbox

    implicit none
    integer :: is

    ! calculate the shock process
    call discretize_AR(rho, 0d0, sigma_eps, eta, pi)

    ! simulate a series of shocks (for both models)
    call simulate_AR(pi, is_sim)

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    do is = 1, ns
        l(:, is) = 0.3d0
        c(:, is) = exp(eta(is))*k(:)**alpha*l(:, is)**(1d0-alpha)
        c_new(:, is) = c(:, is)
    enddo

    ! simulate economy with variable g
    variable_g = .true.
    call getPolicy()
    call output()

    ! store variables for welfare analysis
    WE_c(1)      = E_c
    WE_leis(1)   = 1d0-E_l
    WCV_c(1)     = CV_c
    WCV_leis(1)  = CofVar(1d0-l_t)
    WCor_leis(1) = correlation(c_t, 1d0-l_t)

    ! simulate economy with fixed g
    g_level    = E_g
    variable_g = .false.
    call getPolicy()
    call output()

    ! store variables for welfare analysis
    WE_c(2)      = E_c
    WE_leis(2)   = 1d0-E_l
    WCV_c(2)     = CV_c
    WCV_leis(2)  = CofVar(1d0-l_t)
    WCor_leis(2) = correlation(c_t, 1d0-l_t)

    ! perform welfare analysis
    call welfare_analysis()

contains


    ! Calculates policy functions in the RBS model
    subroutine getPolicy()

        use toolbox

        implicit none
        integer :: ik, is, iter

        ! iterate until policy function converges
        do iter = 1, itermax

            ! interpolate coefficients
            call interpolate()

            ! calculate decisions for every gridpoint and shock level
            do ik = 0, NK
                do is = 1, NS

                    ! set starting value and communicate resource level
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

            ! check for convergence
            if(con_lev < sig)then
                write(*,'(/a,i5,2x,f20.10/)')'CONVERGED: ', iter, con_lev
                return
            endif

            c = c_new
        enddo

        write(*,*)'No Convergence'

    end subroutine


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
                        (1d0+(exp(eta(is_p))*alpha*(k(ik)/l(ik, is_p))**(alpha-1d0)-delta))*marg_util
                enddo
                RHS(ik, is) = (beta*RHS(ik, is))**(-gamma)
            enddo

            ! interpolate
            call spline_interp(RHS(:, is), coeff_r(:, is))
        enddo

    end subroutine


    ! Perform the welfare analysis
    subroutine welfare_analysis

        implicit none
        real*8 :: W_lev(2), W_c(2), W_leis(2), W_cor(2), W_tot(2), util(2)
        real*8 :: Del_lev, Del_c, Del_leis, Del_cor, Del_tot, uhat, bar_nu
        integer :: ii

        ! calculate bar_nu
        bar_nu = (1d0-nu)*nu*gamma

        ! calculate the different welfare components
        do ii = 1, 2
            util(ii)   = (WE_c(ii)**nu*WE_leis(ii)**(1d0-nu))**egam/egam
            W_lev(ii)  = 1d0/(nu*egam)
            W_c(ii)    = -1d0/(2d0*nu*gamma)*(nu**2 + bar_nu)*WCV_c(ii)**2
            W_leis(ii) = -1d0/(2d0*nu*gamma)*((1d0-nu)**2 + bar_nu)*WCV_leis(ii)**2
            W_cor(ii)  = -1d0/(nu*gamma)*(nu*(1d0-nu) - bar_nu)*WCor_leis(ii)*WCV_c(ii)*WCV_leis(ii)
            W_tot(ii)  = W_lev(ii) + W_c(ii) + W_leis(ii) + W_cor(ii)
        enddo

        ! calculate welfare changes due to different components
        uhat     = util(2)/util(1)
        Del_lev  = uhat*W_lev(2)-W_lev(1)
        Del_c    = uhat*W_c(2)-W_c(1)
        Del_leis = uhat*W_leis(2)-W_leis(1)
        Del_cor  = uhat*W_cor(2)-W_cor(1)
        Del_tot  = uhat*W_tot(2)-W_tot(1)

        write(*,'(//a/)')'WELFARE ANALYSIS'
        write(*,'(a)')'                     LEVEL       CV(c)       CV(l)    RHO(c,l)       TOTAL'
        write(*,'(a,4f12.5)')'Values 1:     ',util(1), WCV_c(1), WCV_leis(1), WCor_leis(1)
        write(*,'(a,4f12.5)')'Values 2:     ',util(2), WCV_c(2), WCV_leis(2), WCor_leis(2)
        write(*,'(a,5f12.5)')'Aprox. cost:  ',(/Del_lev, Del_c, Del_leis, Del_cor, Del_tot/)*100d0

    end subroutine


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer :: is

        ! interpolate policy functions
        do is = 1, NS
            call spline_interp(c_new(:, is), coeff_c(:, is))
            call spline_interp(l(:, is), coeff_l(:, is))
        enddo

        ! determine the stochastic steady state
        is_t = floor(dble(NS)/2d0)+1
        call simulate_economy()
        k0 = k_t(TT)

        is_t = is_sim

        ! simulate the economy forward
        call simulate_economy()

        ! calculate business cycle statistics
        call bc_statistics()

        ! write output
        write(*,'(a/)')'BUSINESS CYCLE MOMENTS'
        write(*,'(a)')'                         c       g       i       k       l       y'
        write(*,'(a,5f8.2)')'Mean:             ',(/E_c, E_g, E_i, E_k/)/E_y*100d0, E_l*100d0

        write(*,'(a,6f8.2)')'CV (in %):        ',(/CV_c, CV_g, CV_i, CV_k, CV_l, CV_y/)*100d0
        write(*,'(a,6f8.2)')'Corr. with GDP:   ',Cor_c, Cor_g, Cor_i, Cor_k, Cor_l

    end subroutine


    ! simulates the economy for a given series of shocks
    subroutine simulate_economy()

        use toolbox

        implicit none
        integer :: it

        ! calculate the time path of all economic variables
        eta_t(0) = eta(is_t(0))
        k_t(0) = k0
        c_t(0) = spline_eval(k_t(0), coeff_c(:, is_t(0)), k_l, k_u)
        l_t(0) = spline_eval(k_t(0), coeff_l(:, is_t(0)), k_l, k_u)
        y_t(0) = exp(eta_t(0))*k_t(0)**alpha*l_t(0)**(1d0-alpha)
        if(variable_g)then
            g_t(0) = gy*y_t(0)
        else
            g_t(0) = g_level
        endif
        i_t(0) = y_t(0)-c_t(0)-g_t(0)
        r_t(0) = exp(eta_t(0))*alpha*(l_t(0)/k_t(0))**(1d0-alpha) - delta
        w_t(0) = exp(eta_t(0))*(1d0-alpha)*(k_t(0)/l_t(0))**(alpha)
        do it = 1, TT
            eta_t(it) = eta(is_t(it))
            k_t(it) = (1d0-delta)*k_t(it-1) + i_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c(:, is_t(it)), k_l, k_u)
            l_t(it) = spline_eval(k_t(it), coeff_l(:, is_t(it)), k_l, k_u)
            y_t(it) = exp(eta_t(it))*k_t(it)**alpha*l_t(it)**(1d0-alpha)
            if(variable_g)then
                g_t(it) = gy*y_t(it)
            else
                g_t(it) = g_level
            endif
            i_t(it) = y_t(it)-c_t(it)-g_t(it)
            r_t(it) = exp(eta_t(it))*alpha*(l_t(it)/k_t(it))**(1d0-alpha) - delta
            w_t(it) = exp(eta_t(it))*(1d0-alpha)*(k_t(it)/l_t(it))**alpha
        enddo

    end subroutine


    ! calculates business cycle statistics
    subroutine bc_statistics()

        implicit none

        ! calculate expectations
        E_c = sum(c_t)/dble(TT+1)
        E_g = sum(g_t)/dble(TT+1)
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
        CV_g = CofVar(g_t)
        CV_i = CofVar(i_t)
        CV_k = CofVar(k_t)
        CV_l = CofVar(l_t)
        CV_r = CofVar(r_t)
        CV_w = CofVar(w_t)
        CV_y = CofVar(y_t)

        ! calculate correlation coefficients of logs
        Cor_c = correlation(c_t, y_t)
        Cor_g = correlation(g_t, y_t)
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
