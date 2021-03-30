!##############################################################################
! PROGRAM SOLG_TR_PEN
!
! ## Implicit tax rates in the pension system
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog11_08m.f90"

program SOLG_TR_PEN

    use globals

    implicit none

    ! calculate initial equilibrium
    call get_SteadyState()

    ! set reform parameters (adjust accordingly for Figure 11.8)
    !lambda(0:TT) = 0.99d0
    kappa(1:TT) = 0.19d0

    ! calculate transition path without lsra
    lsra_on = .false.
    call get_transition

    ! calculate transition path with lsra
    lsra_on = .true.
    call get_transition

    ! close files
    close(21)
    close(22)

contains


    ! computes the initial steady state of the economy
    subroutine get_SteadyState()

        implicit none
        integer :: iter

        ! initialize remaining variables
        call initialize()

        ! start timer
        call tic()

        ! iterate until value function converges
        do iter = 1, itermax

            ! derive prices
            call prices(0)

            ! get implicit tax rates
            call implicit_taxes(0)

            ! solve the household problem
            call solve_household(1, 0)

            ! calculate the distribution of households over state space
            call get_distribution(0)

            ! aggregate individual decisions
            call aggregation(0)

            ! determine the government parameters
            call government(0)

            write(*,'(i4,6f8.2,f12.5)')iter, HH(0), (/5d0*KK(0), CC(0), II(0)/)/YY(0)*100d0, &
                ((1d0+r(0))**0.2d0-1d0)*100d0, w(0), DIFF(0)/YY(0)*100d0
            if(abs(DIFF(0)/YY(0))*100d0 < sig)then
                call toc
                call output(0)
                return
            endif
        enddo

        call toc
        call output(0)

        write(*,'(a/)')'ATTENTION: NO CONVERGENCE !!!'

    end subroutine


    ! computes the transition path of the economy
    subroutine get_transition()

        implicit none
        integer :: iter, ij, it, itmax, itcheck
        logical :: check

        ! initialize remaining variables
        if(.not. lsra_on)then
            call initialize_trn()
        else
            write(*,'(a)')'ITER    COMP_OLD  EFFICIENCY        DIFF'
        endif

        ! start timer
        call tic()

        ! iterate until value function converges
        do iter = 1, itermax

            ! derive prices
            do it = 1, TT
                call prices(it)
            enddo

            ! get implicit tax rates
            do it = 1, TT
                call implicit_taxes(it)
            enddo

            ! solve the household problem
            do ij = JJ, 2, -1
                call solve_household(ij, 1)
            enddo
            do it = 1, TT
                call solve_household(1, it)
            enddo

            ! calculate the distribution of households over state space
            do it = 1, TT
                call get_distribution(it)
            enddo

            ! calculate lsra transfers if needed
            if(lsra_on)call LSRA

            ! aggregate individual decisions
            do it = 1, TT
                call aggregation(it)
            enddo

            ! determine the government parameters
            do it = 1, TT
                call government(it)
            enddo

            ! get differences on goods markets
            itcheck = 0
            do it = 1, TT
                if(abs(DIFF(it)/YY(it))*100d0 < sig)itcheck = itcheck + 1
            enddo
            itmax = maxloc(abs(DIFF(1:TT)/YY(1:TT)), 1)

            ! check for convergence and write screen output
            if(.not. lsra_on)then

                check = iter > 1 .and. itcheck == TT-1 .and. &
                    abs(DIFF(itmax)/YY(itmax))*100d0 < sig*100d0

                write(*,'(i4,6f8.2,f12.5)')iter, HH(TT), (/5d0*KK(TT), CC(TT), II(TT)/)/YY(TT)*100d0, &
                    ((1d0+r(TT))**0.2d0-1d0)*100d0, w(TT), DIFF(itmax)/YY(itmax)*100d0
            else
                check = iter > 1 .and. itcheck == TT-1 .and. lsra_comp/lsra_all > 0.99999d0 .and. &
                    abs(DIFF(itmax)/YY(itmax))*100d0 < sig*100d0

                write(*,'(i4,3f12.5)')iter, lsra_comp/lsra_all*100d0, &
                    (Lstar**(1d0/(1d0-1d0/gamma))-1d0)*100d0,DIFF(itmax)/YY(itmax)*100d0
            endif

            ! check for convergence
            if(check)then
                call toc
                do it = 1, TT
                    if(.not. lsra_on)call output(it)
                enddo
                call output_summary()
                return
            endif
        enddo

        call toc
        do it = 1, TT
            if(.not. lsra_on)call output(it)
        enddo
        call output_summary()

        write(*,'(a/)')'ATTENTION: NO CONVERGENCE !!!'

    end subroutine


    ! initializes the remaining model parameters and variables
    subroutine initialize

        implicit none
        integer :: ij, ir, ip, is

        write(*,'(/a/)')'INITIAL EQUILIBRIUM'
        write(*,'(a)')'ITER       H     K/Y     C/Y     I/Y       r       w        DIFF'

        ! set up population structure
        do ij = 1, JJ
            pop(ij, 0) = 1d0/(1d0+n_p)**(ij-1)
        enddo
        do ij = 1, JJ
            m(ij, 0) = pop(ij, 0)/pop(1, 0)
        enddo

        ! initialize asset and pension grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)
        call grid_Cons_Grow(ep, ep_l, ep_u, ep_grow)

        ! get initial guess for savings decision and labor supply
        do ij = 1, JJ
            do ir = 0, NR
                do ip = 1, NP
                    do is = 1, NS
                        aplus(ij, :, ir, ip, is, 0) = max(a(:)/2d0, a(1)/2d0)
                        c(ij, :, ir, ip, is, 0) = max(a(:)/2d0, a(1)/2d0)
                    enddo
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
        tauw  = 0.2088d0
        taur  = 0.2088d0
        taup  = 0.12d0
        kappa = 0.50d0
        !lambda = 1d0
        lambda = 0d0
        penp = 0d0
        gy    = 0.19d0
        by    = 0.60d0/5d0

        ! initial guesses for macro variables
        KK = 5.14d0
        LL = 5.34d0
        YY = 8.42d0
        II = (n_p+delta)*KK
        INC = 0.72d0

        GG = gy*YY(0)
        BB = by*YY(0)
        PP = 0d0

        ! open files
        open(21, file='output.out')
        open(22, file='summary.out')

    end subroutine


    ! initializes transitional variables
    subroutine initialize_trn

        implicit none
        integer :: ij, it

        write(*,'(/a/)')'TRANSITION PATH'

        write(*,'(a)')'ITER       H     K/Y     C/Y     I/Y       r       w        DIFF'

        ! set up population structure
        do it = 1, TT
            pop(1, it) = (1d0+n_p)*pop(1, it-1)
            do ij = 2, JJ
                pop(ij, it) = pop(ij-1, it-1)
            enddo
        enddo

        do it = 1, TT
            do ij = 1, JJ
                m(ij, it) = pop(ij, it)/pop(1, it)
            enddo
        enddo

        do it = 1, TT

            taup(it) = taup(0)
            if(tax(it) == 1)then
                tauc(it) = tauc(0)
            elseif(tax(it) == 2)then
                tauw(it) = tauw(0)
                taur(it) = taur(0)
            elseif(tax(it) == 3)then
                tauw(it) = tauw(0)
            else
                taur(it) = taur(0)
            endif

            r(it) = r(0)
            rn(it) = r(it)*(1d0-taur(it))
            w(it) = w(0)
            wn(it) = w(it)*(1d0-tauw(it)-taup(it))
            p(it) = 1d0 + tauc(it)
            KK(it) = KK(0)
            AA(it) = AA(0)
            BB(it) = BB(0)
            LL(it) = LL(0)
            HH(it) = HH(0)
            YY(it) = YY(0)
            CC(it) = CC(0)
            II(it) = II(0)
            GG(it) = GG(0)
            INC(it) = INC(0)
            pen(:,it) = pen(:, 0)
            penp(:,it,:) = penp(:, 0,:)
            taup(it) = taup(0)
            PP(it) = PP(0)
            taxrev(:,it) = taxrev(:, 0)
            c_coh(:, it) = c_coh(:, 0)
            l_coh(:, it) = l_coh(:, 0)
            y_coh(:, it) = y_coh(:, 0)
            a_coh(:, it) = a_coh(:, 0)
            aplus(:, :, :, :, :, it) = aplus(:, :, :, :, :, 0)
            epplus(:, :, :, :, :, it) = epplus(:, :, :, :, :, 0)
            c(:, :, :, :, :, it) = c(:, :, :, :, :, 0)
            l(:, :, :, :, :, it) = l(:, :, :, :, :, 0)
            phi(:, :, :, :, :, it) = phi(:, :, :, :, :, 0)
            VV(:, :, :, :, :, it) = VV(:, :, :, :, :, 0)
            RHS(:, :, :, :, :, it) = RHS(:, :, :, :, :, 0)
        enddo

    end subroutine


    ! subroutine for calculating prices
    subroutine prices(it)

        implicit none
        integer, intent(in) :: it

        r(it) = Omega*alpha*(KK(it)/LL(it))**(alpha-1d0)-delta
        w(it) = Omega*(1d0-alpha)*(KK(it)/LL(it))**alpha
        rn(it) = r(it)*(1d0-taur(it))
        wn(it) = w(it)*(1d0-tauw(it)-taup(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household(ij_in, it_in)

        implicit none
        integer, intent(in) :: it_in, ij_in
        integer :: ij, ia, ir, ip, ip_max, is, is_max, it
        real*8 :: x_in, wagen, wage, available, v_ind, temp
        logical :: check

        ! get decision in the last period of life
        it = year(it_in, ij_in, JJ)
        do ia = 0, NA
            aplus(JJ, ia, :, :, :, it) = 0d0
            epplus(JJ, ia, :, :, :, it) = 0d0
            l(JJ, ia, :, :, :, it) = 0d0
            do ir = 0, NR
                penp(JJ, it, ir) = kappa(it)*INC(it)*ep(ir)
                c(JJ, ia, ir, :, :, it) = ((1d0+rn(it))*a(ia) + penp(JJ, it, ir) + v(JJ, ia, ir, :, :, it))/p(it)
                VV(JJ, ia, ir, :, :, it) = valuefunc(0d0, 0d0, c(JJ, ia, ir, 1, 1, it), &
                l(JJ, ia, ir, 1, 1, it), JJ, 1, 1, it)
            enddo
        enddo

        ! interpolate individual expectations
        call interpolate(JJ, it)

        do ij = JJ-1, ij_in, -1

            it = year(it_in, ij_in, ij)

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
                if(ij >= JR .and. ia == 0 .and. kappa(it) <= 1d-10)then
                    aplus(ij, ia, :, :, :, it) = 0d0
                    epplus(ij, ia, :, :, :, it) = 0d0
                    c(ij, ia, :, :, :, it) = 0d0
                    l(ij, ia, :, :, :, it) = 0d0
                    VV(ij, ia, :, :, :, it) = valuefunc(0d0, 0d0, 0d0, 0d0, ij, 1, 1, it)
                    cycle
                endif

                do ir = 0, NR
                    do ip = 1, ip_max
                        do is = 1, is_max

                            ! pension system
                            if(ij >= JR)then
                                penp(ij, it, ir) = kappa(it)*INC(it)*ep(ir)
                            else
                                penp(ij, it, ir) = 0d0
                            endif

                            ! get initial guess for the individual choices
                            x_in = aplus(ij, ia, ir, ip, is, it)

                            ! set up communication variables
                            ij_com = ij
                            ia_com = ia
                            ir_com = ir
                            ip_com = ip
                            is_com = is
                            it_com = it

                            ! solve the household problem using rootfinding
                            call fzero(x_in, foc, check)

                            ! write screen output in case of a problem
                            if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ia, ir, ip, is, it

                            ! check for borrowing constraint
                            if(x_in < 0d0)then
                                x_in = 0d0
                                temp = foc(x_in)
                            endif

                            ! copy decisions
                            aplus(ij, ia, ir, ip, is, it) = x_in
                            epplus(ij, ia, ir, ip, is, it) = epplus_com
                            c(ij, ia, ir, ip, is, it) = cons_com
                            l(ij, ia, ir, ip, is, it) = lab_com
                            VV(ij, ia, ir, ip, is, it) = valuefunc(x_in, epplus_com, cons_com, lab_com, ij, ip, is, it)
                        enddo

                        ! copy decision in retirement age
                        if(ij >= JR)then
                            aplus(ij, ia, ir, :, :, it) = aplus(ij, ia, ir, 1, 1, it)
                            epplus(ij, ia, ir, :, :, it) = epplus(ij, ia, ir, 1, 1, it)
                            c(ij, ia, ir, :, :, it) = c(ij, ia, ir, 1, 1, it)
                            l(ij, ia, ir, :, :, it) = l(ij, ia, ir, 1, 1, it)
                            VV(ij, ia, ir, :, :, it) = VV(ij, ia, ir, 1, 1, it)
                        endif
                    enddo
                enddo
            enddo

            ! interpolate individual expectations
            call interpolate(ij, it)
        enddo

    end subroutine


    ! for interpolating individual expectations
    subroutine interpolate(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: ia, ir, ip, is, is_p
        real*8 :: chelp, lhelp

        do ia = 0, NA
            do ir = 0, NR
                do ip = 1, NP
                    do is = 1, NS

                        ! calculate RHS and EV
                        RHS(ij, ia, ir, ip, is, it) = 0d0
                        EV(ij, ia, ir, ip, is, it) = 0d0
                        do is_p = 1, NS
                            chelp = max(c(ij, ia, ir, ip, is_p, it),1d-10)
                            lhelp = max(l(ij, ia, ir, ip, is_p, it),1d-10)
                            RHS(ij, ia, ir, ip, is, it) = RHS(ij, ia, ir, ip, is, it) + &
                                                              pi(is, is_p)*margu(chelp, lhelp, it)
                            EV(ij, ia, ir, ip, is, it)  = EV(ij, ia, ir, ip, is, it) + &
                                                             pi(is, is_p)*VV(ij, ia, ir, ip, is_p, it)
                        enddo
                        RHS(ij, ia, ir, ip, is, it) = ((1d0+rn(it))*beta*RHS(ij, ia, ir, ip, is, it))**(-gamma)
                        EV(ij, ia, ir, ip, is, it) = (egam*EV(ij, ia, ir, ip, is, it))**(1d0/egam)
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, ia, ir, ip, is, is_p, ial, iar, irl, irr, itm
        real*8 :: varphi_a, varphi_r

        ! get yesterdays year
        itm = year(it, 2, 1)

        ! set distribution to zero
        phi(:, :, :, :, :, it) = 0d0

        ! get initial distribution in age 1
        do ip = 1, NP
            phi(1, 0, 0, ip, is_initial, it) = dist_theta(ip)
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS

                            ! interpolate yesterday's savings decision
                            call linint_Grow(aplus(ij-1, ia, ir, ip, is, itm), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                            call linint_Grow(epplus(ij-1, ia, ir, ip, is, itm), ep_l, ep_u, &
                                                                                      ep_grow, NR, irl, irr, varphi_r)

                            ! restrict values to grid just in case
                            ial = min(ial, NA)
                            iar = min(iar, NA)
                            varphi_a = max(0d0, min(varphi_a, 1d0))
                            irl = min(irl, NR)
                            irr = min(irr, NR)
                            varphi_r = max(0d0, min(varphi_r, 1d0))

                            ! redistribute households
                            do is_p = 1, NS
                                phi(ij, ial, irl, ip, is_p, it) = phi(ij, ial, irl, ip, is_p, it) + &
                                                pi(is, is_p)*varphi_a*varphi_r*phi(ij-1, ia, ir, ip, is, itm)
                                phi(ij, ial, irr, ip, is_p, it) = phi(ij, ial, irr, ip, is_p, it) + &
                                                pi(is, is_p)*varphi_a*(1d0-varphi_r)*phi(ij-1, ia, ir, ip, is, itm)
                                phi(ij, iar, irl, ip, is_p, it) = phi(ij, iar, irl, ip, is_p, it) + &
                                                pi(is, is_p)*(1d0-varphi_a)*varphi_r*phi(ij-1, ia, ir, ip, is, itm)
                                phi(ij, iar, irr, ip, is_p, it) = phi(ij, iar, irr, ip, is_p, it) + &
                                          pi(is, is_p)*(1d0-varphi_a)*(1d0-varphi_r)*phi(ij-1, ia, ir, ip, is, itm)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain
    subroutine aggregation(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, ia, ir, ip, is, itp
        real*8 :: LL_old, m_coh(JJ)

        ! get tomorrow's year
        itp = year(it, 1, 2)
        LL_old = LL(it)

        ! calculate cohort aggregates
        c_coh(:, it)  = 0d0
        l_coh(:, it)  = 0d0
        y_coh(:, it)  = 0d0
        a_coh(:, it)  = 0d0
        pen(:, it)    = 0d0
        VV_coh(:, it) = 0d0
        m_coh(:)      = 0d0
        FLC(:,it)     = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS
                            c_coh(ij, it) = c_coh(ij, it) + c(ij, ia, ir, ip, is, it)*phi(ij, ia, ir, ip, is, it)
                            l_coh(ij, it) = l_coh(ij, it) + l(ij, ia, ir, ip, is, it)*phi(ij, ia, ir, ip, is, it)
                            y_coh(ij, it) = y_coh(ij, it) + eff(ij)*theta(ip)*eta(is)*l(ij, ia, ir, ip, is, it) &
                                                                                    *phi(ij, ia, ir, ip, is, it)
                            a_coh(ij, it) = a_coh(ij, it) + a(ia)*phi(ij, ia, ir, ip, is, it)
                            pen(ij, it) = pen(ij, it) +  penp(ij, it, ir)*phi(ij, ia, ir, ip, is, it)
                            ! exclude households who dies
                            if(ij >= JR .and. ia == 0 .and. (kappa(0) <= 1d-10 .or. kappa(1) <= 1d-10))then
                                cycle
                            endif
                            if(aplus(ij, ia, ir, ip, is, it) < 1d-4)FLC(ij, it) = FLC(ij, it) + &
                                                                                    phi(ij, ia, ir, ip, is, it)
                            VV_coh(ij, it) = VV_coh(ij, it) + VV(ij, ia, ir, ip, is, it)*phi(ij, ia, ir, ip, is, it)
                            m_coh(ij)      = m_coh(ij) + phi(ij, ia, ir, ip, is, it)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! normalize VV_coh (because hh excluded)
        VV_coh(:, it) = VV_coh(:, it)/m_coh(:)
        FLC(:, it) = FLC(:, it)/m_coh(:)

        ! calculate aggregate quantities
        CC(it) = 0d0
        LL(it) = 0d0
        HH(it) = 0d0
        AA(it) = 0d0
        workpop(it) = 0d0
        do ij = 1, JJ
            CC(it) = CC(it) + c_coh(ij, it)*m(ij, it)
            LL(it) = LL(it) + y_coh(ij, it)*m(ij, it)
            HH(it) = HH(it) + l_coh(ij, it)*m(ij, it)
            AA(it) = AA(it) + a_coh(ij, it)*m(ij, it)
            if(ij < JR)workpop(it) = workpop(it) + m(ij, it)
        enddo

        ! damping and other quantities
        KK(it) = damp*(AA(it)-BB(it)-BA(it))+(1d0-damp)*KK(it)
        LL(it) = damp*LL(it) + (1d0-damp)*LL_old
        II(it) = (1d0+n_p)*KK(itp) - (1d0-delta)*KK(it)
        YY(it) = Omega * KK(it)**alpha * LL(it)**(1d0-alpha)

        ! get average income and average working hours
        INC(it) = w(it)*LL(it)/workpop(it)
        HH(it)  = HH(it)/workpop(it)

        ! get difference on goods market
        DIFF(it) = YY(it)-CC(it)-II(it)-GG(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, itm, itp
        real*8 :: expend, EP_total

        ! last year
        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        ! set government quantities and pension payments
        GG(it) = gy*YY(0)
        BB(it) = by*YY(0)
        PP(it) = 0d0
        do ij = JR, JJ
            PP(it) = PP(it) + pen(ij, it)*m(ij, it)
        enddo

        ! calculate government expenditure
        expend = GG(it) + (1d0+r(it))*BB(it) - (1d0+n_p)*BB(itp)

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = (expend + - (tauw(it)*w(it)*LL(it) + taur(it)*r(it)*AA(it)))/CC(it)
            p(it)    = 1d0 + tauc(it)
        elseif(tax(it) == 2)then
            tauw(it) = (expend + - tauc(it)*CC(it))/(w(it)*LL(it) + r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = (expend + - (tauc(it)*CC(it) + taur(it)*r(it)*AA(it)))/(w(it)*LL(it))
        else
            taur(it) = (expend + - (tauc(it)*CC(it) + tauw(it)*w(it)*LL(it)))/(r(it)*AA(it))
        endif

        taxrev(1, it) = tauc(it)*CC(it)
        taxrev(2, it) = tauw(it)*w(it)*LL(it)
        taxrev(3, it) = taur(it)*r(it)*AA(it)
        taxrev(4, it) = sum(taxrev(1:3, it))

        ! get budget balancing social security replacement rate
        EP_total = PP(it)/kappa(it)
        taup(it) = kappa(it)*EP_total/(w(it)*LL(it))

    end subroutine


    ! subroutine for calculating implicit taxes in the pension system
    subroutine implicit_taxes(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, ijp, itp

        do ij = 1, JR-1
            itp = year(it, ij, JJ)
            tau_impl(ij, it) = 0d0
            do ijp = JJ, ij+1, -1
                itp = year(it, ij, ijp)
                if(ijp >= JR)then
                    tau_impl(ij, it) = tau_impl(ij, it) + kappa(itp)*INC(itp)
                endif
                tau_impl(ij, it) = tau_impl(ij, it)/(1d0+rn(itp))
            enddo
            tau_impl(ij, it) = taup(it) - (1d0-lambda(it))*tau_impl(ij, it)/(dble(JR-1)*INC(it))
        enddo


    end subroutine


    ! subroutine for calculating lsra payments
    subroutine LSRA()

        implicit none
        integer :: ij, ia, ir, ip, is, it
        real*8 :: VV_1, VV_0, dVV_dv, v_tilde
        real*8 :: EVV_t, EVV_0, dEVV_dv
        real*8 :: PV_t, PV_0, PV_trans

        ! initialize variables
        SV(:) = 0d0
        v_coh(:, :) = 0d0

        ! initialize counters
        lsra_comp     = 0d0
        lsra_all      = 0d0

        do ij = 2, JJ
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS

                            ! do not do anything for an agent at retirement without pension and savings
                            if(ij >= JR .and. ia == 0 .and. (kappa(0) <= 1d-10 .or. kappa(1) <= 1d-10))then
                                v(ij, ia, ir, ip, is, 1) = 0d0
                                cycle
                            endif

                            ! get today's utility
                            VV_1 = VV(ij, ia, ir, ip, is, 1)

                            ! get target utility
                            VV_0 = VV(ij, ia, ir, ip, is, 0)


                            ! get derivative of the value function
                            dVV_dv = margu(c(ij, ia, ir, ip, is, 1),l(ij, ia, ir, ip, is, 1), 1)


                            ! calculate change in transfers
                            v_tilde = (VV_0-VV_1)/dVV_dv

                            ! restrict z_tilde to income maximum
                            v_tilde = max(v_tilde, -((1d0+rn(1))*a(ia) + penp(ij, 1, ir) &
                                      + wn(1)*eff(ij)*theta(ip)*eta(is)*0.99d0 + v(ij, ia, ir, ip, is, 1)))

                            ! check whether individual is already compensated
                            lsra_all = lsra_all + phi(ij, ia, ir, ip, is, 1)*m(ij, 1)
                            if(abs((VV_1-VV_0)/VV_0)*100d0 < sig) &
                                lsra_comp = lsra_comp + phi(ij, ia, ir, ip, is, 1)*m(ij, 1)

                            ! calculate total transfer
                            v(ij, ia, ir, ip, is, 1) = v(ij, ia, ir, ip, is, 1) + damp*v_tilde

                            ! aggregate transfers by cohort
                            v_coh(ij, 1) = v_coh(ij, 1) + v(ij, ia, ir, ip, is, 1)*phi(ij, ia, ir, ip, is, 1)

                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! aggregate transfers in year 1
        do ij = 2, JJ
            SV(1) = SV(1) + v_coh(ij, 1)*m(ij, 1)
        enddo

        ! initialize present value variables
        PV_t = 0d0
        PV_0 = 0d0
        PV_trans = 0d0

        ! calculate present value of utility changes (in monetary values)
        do it = TT, 1, -1

            ! get today's ex ante utility
            EVV_t = damp*VV_coh(1, it)

            ! get damped target utility
            EVV_0 = damp*VV_coh(1, 0)

            ! get derivative of expected utility function
            dEVV_dv = 0d0
            do ip = 1, NP
                do is = 1, NS
                    dEVV_dv = dEVV_dv + margu(c(1, 0, 0, ip, is, it), l(1, 0, 0, ip, is, it), it) &
                                        *phi(1, 0, 0, ip, is, it)
                enddo
            enddo

            ! calculate present values
            if(it == TT)then
                PV_t     = EVV_t/dEVV_dv    *(1d0+r(it))/(r(it)-n_p)
                PV_0     = EVV_0/dEVV_dv    *(1d0+r(it))/(r(it)-n_p)
                PV_trans = v(1, 0, 0, 1, 1, it)*(1d0+r(it))/(r(it)-n_p)
            else
                PV_t     = PV_t    *(1d0+n_p)/(1d0+r(it+1)) + EVV_t/dEVV_dv
                PV_0     = PV_0    *(1d0+n_p)/(1d0+r(it+1)) + EVV_0/dEVV_dv
                PV_trans = PV_trans*(1d0+n_p)/(1d0+r(it+1)) + v(1, 0, 0, 1, 1, it)
            endif
        enddo


        ! calculate the constant utility gain/loss for future generations
        Lstar = (PV_t-PV_trans-SV(1))/PV_0

        ! calculate compensation payments for future cohorts
        do it = TT, 1, -1

            ! get today's ex ante utility
            EVV_t = damp*VV_coh(1, it)

            ! get target utility
            EVV_0 = damp*VV_coh(1, 0)*Lstar

            ! get derivative of expected utility function
            dEVV_dv = 0d0
            do ip = 1, NP
                do is = 1, NS
                    dEVV_dv = dEVV_dv + margu(c(1, 0, 0, ip, is, it), l(1, 0, 0, ip, is, it), it) &
                                        *phi(1, 0, 0, ip, is, it)
                enddo
            enddo

            ! compute change in transfers (restricted)
            v_tilde = (EVV_0-EVV_t)/dEVV_dv

            ! calculate cohort transfer level
            v(1, 0, 0, :, :, it) = v(1, 0, 0, :, :, it) + v_tilde

            ! aggregate transfers
            v_coh(1, it) = v(1, 0, 0, 1, 1, it)
            SV(it) = SV(it) + v_coh(1, it)*m(1, it)

        enddo

        ! determine sequence of LSRA debt/savings
        BA(2) = SV(1)/(1d0+n_p)
        do it = 3, TT
            BA(it) = ((1d0+r(it-1))*BA(it-1) + SV(it-1))/(1d0+n_p)
        enddo

    end subroutine


    ! subroutine for writing output
    subroutine output(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, ia, ir, ip, is, itp, iamax(JJ), iemax(JJ)
        real*8 :: temp, diff_coh
        real*8 :: exp_c(JJ), exp_l(JJ), exp_y(JJ)
        real*8 :: var_c(JJ), var_l(JJ), var_y(JJ)
        real*8 :: mas_c(JJ), mas_l(JJ), mas_y(JJ)

        ! calculate cohort specific variances of logs
        exp_c = 0d0 ; var_c = 0d0 ; mas_c = 0d0
        exp_l = 0d0 ; var_l = 0d0 ; mas_l = 0d0
        exp_y = 0d0 ; var_y = 0d0 ; mas_y = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS

                            ! consumption
                            if(c(ij, ia, ir, ip, is, it) > 0d0)then
                                temp = log(c(ij, ia, ir, ip, is, it))
                                exp_c(ij) = exp_c(ij) + temp*phi(ij, ia, ir, ip, is, it)
                                var_c(ij) = var_c(ij) + temp**2*phi(ij, ia, ir, ip, is, it)
                                mas_c(ij) = mas_c(ij) + phi(ij, ia, ir, ip, is, it)
                            endif

                            if(l(ij, ia, ir, ip, is, it) > 0.01d0)then

                                ! hours
                                temp = log(l(ij, ia, ir, ip, is, it))
                                exp_l(ij) = exp_l(ij) + temp*phi(ij, ia, ir, ip, is, it)
                                var_l(ij) = var_l(ij) + temp**2*phi(ij, ia, ir, ip, is, it)
                                mas_l(ij) = mas_l(ij) + phi(ij, ia, ir, ip, is, it)

                                ! earnings
                                temp = log(w(it)*eff(ij)*theta(ip)*eta(is)*l(ij, ia, ir, ip, is, it))
                                exp_y(ij) = exp_y(ij) + temp*phi(ij, ia, ir, ip, is, it)
                                var_y(ij) = var_y(ij) + temp**2*phi(ij, ia, ir, ip, is, it)
                                mas_y(ij) = mas_y(ij) + phi(ij, ia, ir, ip, is, it)
                            endif
                        enddo
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


        ! Output
        write(21,'(a, i3/)')'EQUILIBRIUM YEAR ', it
        write(21,'(a)')'CAPITAL        K       A       B      BA       r    p.a.'
        write(21,'(8x,6f8.2)')KK(it), AA(it), BB(it), BA(it),  r(it) , ((1d0+r(it))**(1d0/5d0)-1d0)*100d0
        write(21,'(a,4f8.2/)')'(in %)  ',(/KK(it), AA(it), BB(it), BA(it)/)/YY(it)*500d0

        write(21,'(a)')'LABOR          L      HH     INC       w'
        write(21,'(8x,4f8.2/)')LL(it), HH(it)*100d0, INC(it), w(it)

        write(21,'(a)')'GOODS          Y       C       I       G    DIFF'
        write(21,'(8x,4f8.2,f8.3)')YY(it),CC(it),II(it),GG(it),diff(it)
        write(21,'(a,4f8.2,f8.3/)')'(in %)  ',(/YY(it), CC(it), II(it), GG(it), diff(it)/)/YY(it)*100d0

        write(21,'(a)')'GOV         TAUC    TAUW    TAUR    TOTAL      G       B'
        write(21,'(8x,7f8.2)')taxrev(1:4, it),GG(it),BB(it)
        write(21,'(a,6f8.2)')'(in %)  ',(/taxrev(1:4, it), GG(it), BB(it)*5d0/)/YY(it)*100d0
        write(21,'(a,3f8.2/)')'(rate)  ',(/tauc(it), tauw(it), taur(it)/)*100d0

        write(21,'(a)')'PENS        TAUP     PEN      PP'
        write(21,'(8x,6f8.2)')taup(it)*w(it)*LL(it), pen(JR, it), PP(it)
        write(21,'(a,3f8.2/)')'(in %)  ',(/taup(it), kappa(it), PP(it)/YY(it)/)*100d0

        write(21,'(a)')'LSRA          SV      BA'
        write(21,'(8x,2f8.2)')SV(it), BA(it)
        write(21,'(a,2f8.2/)')'(in %)  ',(/SV(it), BA(it)/)/YY(it)*100d0


        ! check for the maximium grid point used
        call check_grid(iamax, iemax, it)
        itp = year(it, 1, 2)

        write(21, '(a,a)')' IJ      CONS     LABOR  EARNINGS    INCOME    INCTAX      PENS    ASSETS', &
            '    VAR(C)    VAR(L)    VAR(Y)      LSRA     VALUE     FLC       IAMAX     IEMAX'
        do ij = 1, JJ
            write(21,'(i3,13f10.3,2i10)')ij, c_coh(ij,it)/INC(0), l_coh(ij,it), (/w(it)*y_coh(ij,it), &
                    wn(it)*y_coh(ij,it)+rn(it)*a_coh(ij,it), &
                    tauw(it)*w(it)*y_coh(ij,it)+taur(it)*r(it)*a_coh(ij,it), &
                    pen(ij,it)-taup(it)*w(it)*y_coh(ij,it), 5d0*a_coh(ij,it)/)/INC(0), &
                    var_c(ij), var_l(ij), var_y(ij), v_coh(ij, it), VV_coh(ij,it), FLC(ij, it), iamax(ij), iemax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

    end subroutine


    ! writes summary output
    subroutine output_summary()

        implicit none
        integer :: ij, ia, ir, ip, is, it
        real*8 :: HEV(-(JJ-2):TT), HEV_help, mas(-(JJ-2):0)

        ! aggregate ex post welfare changes of current generations
        HEV = 0d0
        mas = 0d0
        do ij = JJ, 2, -1
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS
                            if(ij >= JR .and. ia == 0 .and. (kappa(0) <= 1d-10 .or. kappa(1) <= 1d-10))then
                                cycle
                            endif
                            HEV_help = ((VV(ij, ia, ir, ip, is, 1) &
                                      /max(VV(ij, ia, ir, ip, is, 0), -1d10))**(1d0/egam)-1d0)*100d0
                            HEV(-(ij-2)) = HEV(-(ij-2)) + HEV_help*phi(ij, ia, ir, ip, is, 1)
                            mas(-(ij-2)) = mas(-(ij-2)) + phi(ij, ia, ir, ip, is, 1)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        HEV(-(JJ-2):0) = HEV(-(JJ-2):0)/mas

        ! calculate ex ante welfare of future generations
        do it = 1, TT
            HEV(it) = ((VV_coh(1, it)/VV_coh(1, 0))**(1d0/egam)-1d0)*100d0
        enddo

        ! headline
        write(22,'(/a,a)')'          A       K       L       H       r       w       C       I       Y', &
            '       B      BA    tauc    tauw    taur    taup     HEV      DIFF'

        ! current generations
        do ij = -(JJ-2), -1
            write(22,'(i3,120x,f8.2)')ij,HEV(ij)
        enddo

        ! future generations
        do it = 0, TT
            write(22,'(i3,16f8.2,f10.5)')it,(/AA(it)/AA(0)-1d0, KK(it)/KK(0)-1d0, LL(it)/LL(0)-1d0, &
                HH(it)-HH(0), (1d0+r(it))**0.2d0-(1d0+r(0))**0.2d0, w(it)/w(0)-1d0, CC(it)/CC(0)-1d0, &
                II(it)/II(0)-1d0, YY(it)/YY(0)-1d0, BB(it)/BB(0)-1d0, BA(it)/YY(it), tauc(it)-tauc(0), &
                tauw(it)-tauw(0), taur(it)-taur(0), taup(it)-taup(0)/)*100d0, &
                HEV(it), DIFF(it)/YY(it)*100d0
        enddo

        if(lsra_on)write(22, '(/a,f12.6)')'EFFICIENCY GAIN: ', (Lstar**(1d0/(1d0-1d0/gamma))-1d0)*100d0

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax, iemax, it)

        implicit none
        integer, intent(in) :: it
        integer :: iamax(JJ), iemax(JJ), ij, ia, ir, ip, is

        iamax = 0
        iemax = 0

        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do ir = 0, NR
                    do ip = 1, NP
                        do is = 1, NS
                            if(phi(ij, ia, ir, ip, is, it) > 1d-6)iamax(ij) = ia
                        enddo
                    enddo
                enddo
            enddo

            ! check for the maximum earning point grid point used at a certain age
            do ir = 0, NR
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS
                            if(phi(ij, ia, ir, ip, is, it) > 1d-6)iemax(ij) = ir
                        enddo
                    enddo
                enddo
            enddo

        enddo

    end subroutine

end program
