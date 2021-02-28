!##############################################################################
! PROGRAM SOLG_LR
!
! ## Long-run equilibria in the stochastic OLG model
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
include "prog11_01m.f90"

program SOLG_LR

    use globals

    implicit none

    ! initialize variables
    call initialize()

    open(21, file='output_initial.out')

    ! calculate initial equilibrium
    reform_on = .false.
    call get_SteadyState()

    close(21)

    open(21, file='output_final.out')

    ! set reform variables
    reform_on = .true.

    ! set reform values (Table 11.3, different rows)
    tax = 1    ;    tauw = 0d0    ;    taur = 0d0  ! Simulation (1)
    !tax = 3    ;    taur = 0d0                     ! Simulation (2)
    !kappa = 0d0                                    ! Simulation (3)

    ! calculate final equilibrium
    call get_SteadyState()

    close(21)

contains


    ! computes the initial steady state of the economy
    subroutine get_SteadyState()

        implicit none
        integer :: iter

        ! start timer
        call tic()

        ! iterate until value function converges
        do iter = 1, itermax

            ! derive prices
            call prices()

            ! solve the household problem
            call solve_household()

            ! calculate the distribution of households over state space
            call get_distribution()

            ! aggregate individual decisions over cohorts
            call aggregation()

            ! determine the government parameters
            call government()

            write(*,'(i4,5f8.2,f12.5)')iter, (/5d0*KK, CC, II/)/YY*100d0, &
                                       r, w, DIFF/YY*100d0
            if(abs(DIFF/YY)*100d0 < sig)then
                call toc
                call output()
                return
            endif
        enddo

        call toc
        call output()

        write(*,'(a/)')'ATTENTION: NO CONVERGENCE !!!'

    end subroutine


    ! initializes the remaining model parameters and variables
    subroutine initialize

        implicit none
        integer :: ij, ip, is

        write(*,'(/a/)')'INITIAL EQUILIBRIUM'
        write(*,'(a)')'ITER     K/Y     C/Y     I/Y       r       w        DIFF'

        ! set up population structure
        do ij = 1, JJ
            m(ij) = (1d0+n_p)**(1d0-ij)
        enddo

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! get initial guess for savings decision
        do ij = 1, JJ
            do ip = 1, NP
                do is = 1, NS
                    aplus(ij, :, ip, is) = max(a(:)/2d0, a(1)/2d0)
                enddo
            enddo
        enddo

        ! initialize age earnings process
        eff(1) = 1.0000d0
        eff(2) = 1.3527d0
        eff(3) = 1.6952d0
        eff(4) = 1.8279d0
        eff(5) = 1.9606d0
        eff(6) = 1.9692d0
        eff(7) = 1.9692d0
        eff(8) = 1.9392d0
        eff(9) = 1.9007d0
        eff(JR:JJ) = 0d0

        ! initialize fixed effect
        dist_theta = 1d0/dble(NP)
        theta(1)   = -sqrt(sigma_theta)
        theta(2)   = sqrt(sigma_theta)
        theta = exp(theta)

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! tax and transfers
        tax   = 2
        tauc  = 0.075d0
        tauw  = 0.0d0
        taur  = 0.0d0
        taup  = 0.1d0
        kappa = 0.5d0
        gy    = 0.19d0
        by    = 0.60d0/5d0

        ! initial guesses for macro variables
        KK = 1d0
        LL = 1d0
        YY = 1d0
        II = (n_p+delta)*KK

        GG = gy*YY
        BB = by*YY

        pen = 0d0
        pen(JR:JJ) = kappa

    end subroutine


    ! subroutine for calculating prices
    subroutine prices()

        implicit none

        r = Omega*alpha*(KK/LL)**(alpha-1d0)-delta
        w = Omega*(1d0-alpha)*(KK/LL)**alpha
        rn = r*(1d0-taur)
        wn = w*(1d0-tauw-taup)
        p = 1d0 + tauc

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, ip, ip_max, is, is_max
        real*8 :: x_in, wage, available
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            aplus(JJ, ia, :, :) = 0d0
            c(JJ, ia, :, :) = ((1d0+rn)*a(ia) + pen(JJ))/p
            l(JJ, ia, :, :) = 0d0
            V(JJ, ia, :, :) = valuefunc(0d0, c(JJ, ia, 1, 1), l(JJ, ia, 1, 1), JJ, 1, 1)
        enddo

        ! interpolate individual RHS
        call interpolate(JJ)

        do ij = JJ-1, 1, -1

            ! check about how many is to iterate
            if(ij >= JR)then
                ip_max = 1
                is_max = 1
            else
                ip_max = NP
                is_max = NS
            endif

            do ia = 0, NA

                ! determine decision for zero assets at retirement without pension
                if(ij >= JR .and. ia == 0 .and. kappa <= 1d-10)then
                    aplus(ij, ia, :, :) = 0d0
                    c(ij, ia, :, :) = 0d0
                    l(ij, ia, :, :) = 0d0
                    V(ij, ia, :, :) = valuefunc(0d0, 0d0, 0d0, ij, 1, 1)
                    cycle
                endif

                do ip = 1, ip_max
                    do is = 1, is_max

                        ! get initial guess for the individual choices
                        x_in = aplus(ij, ia, ip, is)

                        ! set up communication variables
                        ij_com = ij
                        ia_com = ia
                        ip_com = ip
                        is_com = is

                        ! solve the household problem using rootfinding
                        call fzero(x_in, foc, check)

                        ! write screen output in case of a problem
                        if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is

                        ! check for borrowing constraint
                        if(x_in < 0d0)then
                            x_in = 0d0
                            wage = wn*eff(ij)*theta(ip)*eta(is)
                            available = (1d0+rn)*a(ia) + pen(ij)
                            if(ij < JR)then
                                lab_com = min( max(nu-(1d0-nu)*available/wage , 0d0) , 1d0-1d-10)
                            else
                                lab_com = 0d0
                            endif
                            cons_com = max( (available + wage*lab_com)/p , 1d-10)
                        endif

                        ! copy decisions
                        aplus(ij, ia, ip, is) = x_in
                        c(ij, ia, ip, is) = cons_com
                        l(ij, ia, ip, is) = lab_com
                        V(ij, ia, ip, is) = valuefunc(x_in, cons_com, lab_com, ij, ip, is)

                    enddo

                    ! copy decision in retirement age
                    if(ij >= JR)then
                        aplus(ij, ia, :, :) = aplus(ij, ia, 1, 1)
                        c(ij, ia, :, :) = c(ij, ia, 1, 1)
                        l(ij, ia, :, :) = l(ij, ia, 1, 1)
                        V(ij, ia, :, :) = V(ij, ia, 1, 1)
                    endif
                enddo
            enddo

            ! interpolate individual RHS
            call interpolate(ij)
        enddo

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, ip, is, is_p
        real*8 :: chelp, lhelp

        do ia = 0, NA
            do ip = 1, NP
                do is = 1, NS

                    ! calculate the RHS of the first order condition
                    RHS(ij, ia, ip, is) = 0d0
                    EV(ij, ia, ip, is) = 0d0
                    do is_p = 1, NS
                        chelp = max(c(ij, ia, ip, is_p),1d-10)
                        lhelp = max(l(ij, ia, ip, is_p),1d-10)
                        RHS(ij, ia, ip, is) = RHS(ij, ia, ip, is) + pi(is, is_p)*margu(chelp, lhelp)
                        EV(ij, ia, ip, is)  = EV(ij, ia, ip, is) + pi(is, is_p)*V(ij, ia, ip, is_p)
                    enddo
                    RHS(ij, ia, ip, is) = ((1d0+rn)*beta*RHS(ij, ia, ip, is))**(-gamma)
                    EV(ij, ia, ip, is) = (egam*EV(ij, ia, ip, is))**(1d0/egam)
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, ip, is, is_p, ial, iar
        real*8 :: varphi

        ! set distribution to zero
        phi(:, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ip = 1, NP
            phi(1, 0, ip, is_initial) = dist_theta(ip)
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do ip = 1 , NP
                    do is = 1, NS

                        ! interpolate yesterday's savings decision
                        call linint_Grow(aplus(ij-1, ia, ip, is), a_l, a_u, a_grow, NA, ial, iar, varphi)

                        ! restrict values to grid just in case
                        ial = min(ial, NA)
                        iar = min(iar, NA)
                        varphi = min(varphi, 1d0)

                        ! redistribute households
                        do is_p = 1, NS
                            phi(ij, ial, ip, is_p) = phi(ij, ial, ip, is_p) + &
                                                        pi(is, is_p)*varphi*phi(ij-1, ia, ip, is)
                            phi(ij, iar, ip, is_p) = phi(ij, iar, ip, is_p) + &
                                                        pi(is, is_p)*(1d0-varphi)*phi(ij-1, ia, ip, is)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ia, ip, is
        real*8 :: workpop, LL_old

        LL_old = LL

        ! calculate cohort aggregates
        c_coh(:) = 0d0
        l_coh(:) = 0d0
        y_coh(:) = 0d0
        a_coh(:) = 0d0
        v_coh(:) = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        c_coh(ij) = c_coh(ij) + c(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        l_coh(ij) = l_coh(ij) + l(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        y_coh(ij) = y_coh(ij) + eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        a_coh(ij) = a_coh(ij) + a(ia)*phi(ij, ia, ip, is)
                        v_coh(ij) = v_coh(ij) + V(ij, ia, ip, is)*phi(ij, ia, ip, is)
                    enddo
                enddo
            enddo
        enddo

        ! calculate aggregate quantities
        CC = 0d0
        LL = 0d0
        HH = 0d0
        AA = 0d0
        workpop = 0d0
        do ij = 1, JJ
            CC = CC + c_coh(ij)*m(ij)
            LL = LL + y_coh(ij)*m(ij)
            HH = HH + l_coh(ij)*m(ij)
            AA = AA + a_coh(ij)*m(ij)
            if(ij < JR)workpop = workpop + m(ij)
        enddo

        ! damping and other quantities
        KK = damp*(AA-BB) + (1d0-damp)*KK
        LL = damp*LL + (1d0-damp)*LL_old
        II = (n_p+delta)*KK
        YY = Omega * KK**alpha * LL**(1d0-alpha)

        ! get average income and average working hours
        INC = w*LL/workpop
        HH  = HH/workpop

        ! get difference on goods market
        DIFF = YY-CC-II-GG

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government()

        implicit none
        integer :: ij
        real*8 :: expend

        ! set government quantities and pension payments
        if(.not. reform_on)then
            GG = gy*YY
            BB = by*YY
        endif

        ! calculate government expenditure
        expend = GG + (1d0+r)*BB - (1d0+n_p)*BB

        ! get budget balancing tax rate
        if(tax == 1)then
            tauc = (expend - (tauw*w*LL + taur*r*AA))/CC
            p    = 1d0 + tauc
        elseif(tax == 2)then
            tauw = (expend - tauc*CC)/(w*LL + r*AA)
            taur = tauw
        elseif(tax == 3)then
            tauw = (expend - (tauc*CC + taur*r*AA))/(w*LL)
        else
            taur = (expend - (tauc*CC + tauw*w*LL))/(r*AA)
        endif

        taxrev(1) = tauc*CC
        taxrev(2) = tauw*w*LL
        taxrev(3) = taur*r*AA
        taxrev(4) = sum(taxrev(1:3))

        ! get budget balancing social security contribution
        pen(JR:JJ) = kappa*INC
        PP = 0d0
        do ij = JR, JJ
            PP = PP + pen(ij)*m(ij)
        enddo

        taup = PP/(w*LL)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ia, ip, is, iamax(JJ)
        real*8 :: temp
        real*8 :: exp_c(JJ), exp_l(JJ), exp_y(JJ)
        real*8 :: var_c(JJ), var_l(JJ), var_y(JJ)
        real*8 :: mas_c(JJ), mas_l(JJ), mas_y(JJ)

        ! calculate cohort specific variances of logs
        exp_c = 0d0 ; var_c = 0d0 ; mas_c = 0d0
        exp_l = 0d0 ; var_l = 0d0 ; mas_l = 0d0
        exp_y = 0d0 ; var_y = 0d0 ; mas_y = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS

                        ! consumption
                        if(c(ij, ia, ip, is) > 0d0)then
                            temp = log(c(ij, ia, ip, is))
                            exp_c(ij) = exp_c(ij) + temp*phi(ij, ia, ip, is)
                            var_c(ij) = var_c(ij) + temp**2*phi(ij, ia, ip, is)
                            mas_c(ij) = mas_c(ij) + phi(ij, ia, ip, is)
                        endif

                        if(l(ij, ia, ip, is) > 0.01d0)then

                            ! hours
                            temp = log(l(ij, ia, ip, is))
                            exp_l(ij) = exp_l(ij) + temp*phi(ij, ia, ip, is)
                            var_l(ij) = var_l(ij) + temp**2*phi(ij, ia, ip, is)
                            mas_l(ij) = mas_l(ij) + phi(ij, ia, ip, is)

                            ! earnings
                            temp = log(w*eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is))
                            exp_y(ij) = exp_y(ij) + temp*phi(ij, ia, ip, is)
                            var_y(ij) = var_y(ij) + temp**2*phi(ij, ia, ip, is)
                            mas_y(ij) = mas_y(ij) + phi(ij, ia, ip, is)
                        endif
                    enddo
                enddo
            enddo
        enddo
        exp_c = exp_c/max(mas_c, 1d-4) ; var_c = var_c/max(mas_c, 1d-4)
        exp_l = exp_l/max(mas_l, 1d-4) ; var_l = var_l/max(mas_l, 1d-4)
        exp_y = exp_y/max(mas_y, 1d-4) ; var_y = var_y/max(mas_y, 1d-4)
        var_c = var_c - exp_c**2
        var_l = var_l - exp_l**2
        var_y = var_y - exp_y**2

        ! save initial equilibrium average income if no reform
        if(.not. reform_on)INC_init = INC

        ! Output
        write(21,'(a/)')'STEADY STATE EQUILIBRIUM'
        write(21,'(a)')'CAPITAL        K       A       B       r    p.a.'
        write(21,'(8x,5f8.2)')KK, AA, BB,  r , ((1d0+r)**(1d0/5d0)-1d0)*100d0
        write(21,'(a,3f8.2/)')'(in %)  ',(/KK, AA, BB/)/YY*500d0

        write(21,'(a)')'LABOR          L      HH     INC       w'
        write(21,'(8x,4f8.2/)')LL, HH*100d0, INC, w

        write(21,'(a)')'GOODS          Y       C       I       G    DIFF'
        write(21,'(8x,4f8.2,f8.3)')YY,CC,II,GG,diff
        write(21,'(a,4f8.2,f8.3/)')'(in %)  ',(/YY, CC, II, GG, diff/)/YY*100d0

        write(21,'(a)')'GOV         TAUC    TAUW    TAUR   TOTAL       G       B'
        write(21,'(8x,6f8.2)')taxrev(1:4),GG,BB
        write(21,'(a,6f8.2)')'(in %)  ',(/taxrev(1:4), GG, BB*5d0/)/YY*100d0
        write(21,'(a,3f8.2/)')'(rate)  ',(/tauc, tauw, taur/)*100d0

        write(21,'(a)')'PENS        TAUP     PEN      PP'
        write(21,'(8x,6f8.2)')taup*w*LL, pen(JR), PP
        write(21,'(a,3f8.2/)')'(in %)  ',(/taup, kappa, PP/YY/)*100d0

        ! check for the maximium grid point used
        call check_grid(iamax)

        write(21, '(a,a)')' IJ      CONS     LABOR  EARNINGS    INCOME    INCTAX      PENS    ASSETS', &
            '    VAR(C)    VAR(L)    VAR(Y)     VALUE     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,11f10.3,i10)')ij, c_coh(ij)/INC_init, l_coh(ij), (/w*y_coh(ij), wn*y_coh(ij)+rn*a_coh(ij), &
                    tauw*w*y_coh(ij)+taur*r*a_coh(ij), pen(ij)-taup*w*y_coh(ij), 5d0*a_coh(ij)/)/INC_init, &
                    var_c(ij), var_l(ij), var_y(ij), v_coh(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, ip, is

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        if(phi(ij, ia, ip, is) > 1d-8)iamax(ij) = ia
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
