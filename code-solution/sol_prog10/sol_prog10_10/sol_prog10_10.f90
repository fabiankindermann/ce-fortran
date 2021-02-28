!##############################################################################
! PROGRAM LaborSupplySeparableTax
!
! ## Endogenous Labor Supply with progressive taxes
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
include "sol_prog10_10m.f90"

program LaborSupplySeparableTax

    use globals

    implicit none
    integer :: ic

    ! initialize remaining variables
    call initialize()

    do ic = 1, 2

        ! start the clock
        call tic()

        ! solve the household problem with policy function iteration
        if(ic == 1)call solve_household_policy()

        ! solve the household problem with endogenous gridpoints
        if(ic == 2)call solve_household_endogenous()

        ! calculate the distribution of households over state space
        call get_distribution()

        ! aggregate individual decisions
        call aggregation()

        ! stop the clock
        call toc()

        call output()

    enddo

    ! close files
    close(21)

contains


    ! initializes all remaining variables
    subroutine initialize

        implicit none
        integer :: ij

        ! net prices (after taxes and tranfers)
        r = 0.04d0
        w = 1.0d0

        ! set survival probabilities
        psi = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
                0.99906d0, 0.99908d0, 0.99906d0, 0.99907d0, 0.99901d0, &
                0.99899d0, 0.99896d0, 0.99893d0, 0.99890d0, 0.99887d0, &
                0.99886d0, 0.99878d0, 0.99871d0, 0.99862d0, 0.99853d0, &
                0.99841d0, 0.99835d0, 0.99819d0, 0.99801d0, 0.99785d0, &
                0.99757d0, 0.99735d0, 0.99701d0, 0.99676d0, 0.99650d0, &
                0.99614d0, 0.99581d0, 0.99555d0, 0.99503d0, 0.99471d0, &
                0.99435d0, 0.99393d0, 0.99343d0, 0.99294d0, 0.99237d0, &
                0.99190d0, 0.99137d0, 0.99085d0, 0.99000d0, 0.98871d0, &
                0.98871d0, 0.98721d0, 0.98612d0, 0.98462d0, 0.98376d0, &
                0.98226d0, 0.98062d0, 0.97908d0, 0.97682d0, 0.97514d0, &
                0.97250d0, 0.96925d0, 0.96710d0, 0.96330d0, 0.95965d0, &
                0.95619d0, 0.95115d0, 0.94677d0, 0.93987d0, 0.93445d0, &
                0.92717d0, 0.91872d0, 0.91006d0, 0.90036d0, 0.88744d0, &
                0.87539d0, 0.85936d0, 0.84996d0, 0.82889d0, 0.81469d0, &
                0.79705d0, 0.78081d0, 0.76174d0, 0.74195d0, 0.72155d0, &
                0.00000d0/)

        ! initialize age earnings process
        eff(1:JR-1) = &
            (/1.0000d0, 1.0719d0, 1.1438d0, 1.2158d0, 1.2842d0, 1.3527d0, &
              1.4212d0, 1.4897d0, 1.5582d0, 1.6267d0, 1.6952d0, 1.7217d0, &
              1.7438d0, 1.7748d0, 1.8014d0, 1.8279d0, 1.8545d0, 1.8810d0, &
              1.9075d0, 1.9341d0, 1.9606d0, 1.9623d0, 1.9640d0, 1.9658d0, &
              1.9675d0, 1.9692d0, 1.9709d0, 1.9726d0, 1.9743d0, 1.9760d0, &
              1.9777d0, 1.9700d0, 1.9623d0, 1.9546d0, 1.9469d0, 1.9392d0, &
              1.9315d0, 1.9238d0, 1.9161d0, 1.9084d0, 1.9007d0, 1.8354d0, &
              1.7701d0, 1.7048d0/)
        eff(JR:JJ) = 0d0

        ! compute cumulated survival probabilities
        psi_cum(:) = 0d0
        psi_cum(1) = psi(1)
        do ij = 2, JJ+1
            psi_cum(ij) = psi_cum(ij-1)*psi(ij)
        enddo

        ! old-age transfers
        pen = 0d0
        pen(JR:JJ) = 0.5d0*sum(eff)/dble(JR-1)*0.33d0

        ! initialize fixed effect
        dist_theta = 1d0/dble(NP)
        theta(1)   = exp(-sqrt(sigma_theta))
        theta(2)   = exp( sqrt(sigma_theta))

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household_policy()

        implicit none
        integer :: ij, ia, ip, ip_max, is, is_max
        real*8 :: x_in
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            aplus(JJ, ia, :, :) = 0d0
            c(JJ, ia, :, :) = (1d0+r)*a(ia) + pen(JJ)
            l(JJ, ia, :, :) = 0d0
            V_c(JJ, ia, :, :) = valuefunc_c(0d0, c(JJ, ia, 1, 1), JJ, 1, 1)
            V_l(JJ, ia, :, :) = valuefunc_l(0d0, l(JJ, ia, 1, 1), JJ, 1, 1)
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
                if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                    aplus(ij, ia, :, :) = 0d0
                    c(ij, ia, :, :) = 0d0
                    l(ij, ia, :, :) = 0d0
                    V_c(ij, ia, :, :) = valuefunc_c(0d0, 0d0, ij, 1, 1)
                    V_l(ij, ia, :, :) = valuefunc_l(0d0, 0d0, ij, 1, 1)
                    cycle
                endif

                do ip = 1, ip_max
                    do is = 1, is_max

                        ! get initial guess for the individual choices
                        x_in = c(ij+1, ia, ip, is)

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
                        if(a_plus_com < 0d0)then
                            if(ij >= JR)then
                                a_plus_com = a_l
                                lab_com = 0d0
                                x_in = (1d0+r)*a(ia) + pen(ij)
                            else
                                a_plus_com = a_l

                                ! solve the household problem using rootfinding
                                call fzero(x_in, foc_liq, check)

                                ! write screen output in case of a problem
                                if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is
                            endif
                        endif

                        ! copy decisions
                        aplus(ij, ia, ip, is) = a_plus_com
                        c(ij, ia, ip, is) = x_in
                        l(ij, ia, ip, is) = lab_com
                        V_c(ij, ia, ip, is) = valuefunc_c(aplus(ij, ia, ip, is), c(ij, ia, ip, is), ij, ip, is)
                        V_l(ij, ia, ip, is) = valuefunc_l(aplus(ij, ia, ip, is), l(ij, ia, ip, is), ij, ip, is)
                    enddo
                enddo

                ! copy decision in retirement age
                if(ij >= JR)then
                    aplus(ij, ia, :, :) = aplus(ij, ia, 1, 1)
                    c(ij, ia, :, :) = c(ij, ia, 1, 1)
                    l(ij, ia, :, :) = l(ij, ia, 1, 1)
                    V_c(ij, ia, :, :) = V_c(ij, ia, 1, 1)
                    V_l(ij, ia, :, :) = V_l(ij, ia, 1, 1)
                endif
            enddo

            ! interpolate individual RHS
            call interpolate(ij)

            write (*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household_endogenous()

        implicit none
        integer :: ij, ia, ip, ip_max, is, is_max
        real*8 :: wage, x_in, expo
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            aplus(JJ, ia, :, :) = 0d0
            c(JJ, ia, :, :) = (1d0+r)*a(ia) + pen(JJ)
            l(JJ, ia, :, :) = 0d0
            V_c(JJ, ia, :, :) = valuefunc_c(0d0, c(JJ, ia, 1, 1), JJ, 1, 1)
            V_l(JJ, ia, :, :) = valuefunc_l(0d0, l(JJ, ia, 1, 1), JJ, 1, 1)
        enddo

        ! interpolate individual RHS and EVs
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

            do ip = 1, ip_max
                do is = 1, is_max
                    do ia = 0, NA

                        ! get wage
                        wage = w*eff(ij)*theta(ip)*eta(is)

                        ! derive consumption and labor supply
                        c_endo(ia) = RHS(ij+1, ia, ip, is)
                        if(ij >= JR)then
                            l_endo(ia) = 0d0
                        else
                            ! get exponent
                            expo = 1d0/(1d0 + 1d0/chi-varrho)

                            ! get labour supply
                            l_endo(ia) = max(((1d0-tau)*varrho*wage**varrho/nu)**expo*c_endo(ia)**(-expo/gamma), 0d0)
                        endif
                        a_endo(ia) = (c_endo(ia) + a(ia) - pen(ij) - (1d0-tau)*(wage*l_endo(ia))**varrho)/(1d0+r)
                    enddo

                    do ia = 0, NA

                        ! check for liquidity constraint
                        if(a(ia) <= a_endo(0))then

                            if(ij >= JR)then
                                c(ij, ia, ip, is) = (1d0+r)*a(ia) + pen(ij)
                                l(ij, ia, ip, is) = 0d0
                                aplus(ij, ia, ip, is) = 0d0
                            else

                                ! solve first order condition under liquidity constriant
                                ij_com = ij
                                ia_com = ia
                                ip_com = ip
                                is_com = is

                                x_in = c_endo(0)

                                ! solve the household problem using rootfinding
                                call fzero(x_in, foc_liq, check)

                                ! write screen output in case of a problem
                                if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is

                                c(ij, ia, ip, is) = x_in
                                l(ij, ia, ip, is) = lab_com
                                aplus(ij, ia, ip, is) = 0d0
                            endif

                        else

                            c(ij, ia, ip, is) = linint_Gen(a(ia), a_endo, c_endo, ia)
                            if(ij >= JR)then
                                l(ij, ia, ip, is) = 0d0
                            else
                                l(ij, ia, ip, is) = linint_Gen(a(ia), a_endo, l_endo, ia)
                            endif
                            aplus(ij, ia, ip, is) = (1d0+r)*a(ia) + (1d0-tau)*(wage*l(ij, ia, ip, is))**varrho &
                                                    + pen(ij) - c(ij, ia, ip, is)

                        endif

                        ! get value function
                        V_c(ij, ia, ip, is) = valuefunc_c(aplus(ij, ia, ip, is), c(ij, ia, ip, is), ij, ip, is)
                        V_l(ij, ia, ip, is) = valuefunc_l(aplus(ij, ia, ip, is), l(ij, ia, ip, is), ij, ip, is)
                    enddo
                enddo
            enddo

            ! copy decision in retirement age
            if(ij >= JR)then
                do ia = 0, NA
                    aplus(ij, ia, :, :) = aplus(ij, ia, 1, 1)
                    c(ij, ia, :, :) = c(ij, ia, 1, 1)
                    l(ij, ia, :, :) = l(ij, ia, 1, 1)
                    V_c(ij, ia, :, :) = V_c(ij, ia, 1, 1)
                    V_l(ij, ia, :, :) = V_l(ij, ia, 1, 1)
                enddo
            endif

            ! interpolate individual RHS
            call interpolate(ij)

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, ip, is, is_p
        real*8 :: chelp

        do ia = 0, NA
            do ip = 1, NP
                do is = 1, NS

                    ! calculate the RHS of the first order condition
                    RHS(ij, ia, ip, is) = 0d0
                    EV_c(ij, ia, ip, is) = 0d0
                    EV_l(ij, ia, ip, is) = 0d0
                    do is_p = 1, NS
                        chelp = max(c(ij, ia, ip, is_p), 1d-10)
                        RHS(ij, ia, ip, is) = RHS(ij, ia, ip, is) + &
                                              pi(is, is_p)*margu(chelp)
                        EV_c(ij, ia, ip, is) = EV_c(ij, ia, ip, is) + &
                                               pi(is, is_p)*V_c(ij, ia, ip, is_p)
                        EV_l(ij, ia, ip, is) = EV_l(ij, ia, ip, is) + &
                                               pi(is, is_p)*V_l(ij, ia, ip, is_p)
                    enddo
                    RHS(ij, ia, ip, is) = ((1d0+r)*beta*psi(ij)*RHS(ij, ia, ip, is))**(-gamma)
                    EV_c(ij, ia, ip, is) = (egam*EV_c(ij, ia, ip, is))**(1d0/egam)
                    EV_l(ij, ia, ip, is) = (echi*EV_l(ij, ia, ip, is))**(1d0/echi)
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
                do ip = 1, NP
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

        ! calculate cohort averages
        c_coh(:) = 0d0
        n_coh(:) = 0d0
        l_coh(:) = 0d0
        h_coh(:) = 0d0
        a_coh(:) = 0d0
        v_coh(:) = 0d0
        t_coh(:) = 0d0
        tax_revenue = 0d0
        labor_supply = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        c_coh(ij) = c_coh(ij) + c(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        n_coh(ij) = n_coh(ij) + l(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        l_coh(ij) = l_coh(ij) + eff(ij)*theta(ip)*eta(is)* &
                                    l(ij, ia, ip, is)*phi(ij, ia, ip, is)
                        t_coh(ij) = t_coh(ij) + (eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is) - &
                                                 (1d0-tau)*(eff(ij)*theta(ip)*eta(is)*l(ij, ia, ip, is))**varrho) &
                                    *phi(ij, ia, ip, is)
                        h_coh(ij) = h_coh(ij) + eff(ij)*theta(ip)*eta(is)*phi(ij, ia, ip, is)
                        a_coh(ij) = a_coh(ij) + a(ia)*phi(ij, ia, ip, is)
                        v_coh(ij) = v_coh(ij) + (V_c(ij, ia, ip, is) - V_l(ij, ia, ip, is))*phi(ij, ia, ip, is)
                    enddo
                enddo
            enddo
            if(ij < JR)then
                tax_revenue = tax_revenue + psi(ij)*t_coh(ij)
                labor_supply = labor_supply + psi(ij)*l_coh(ij)
            endif
        enddo

        ! print tax revenue and labor supply to the screen
        write (*,'(a,2f10.2)')'Tax revenue, Labor supply =', tax_revenue, labor_supply

        ! calculate cohort specific coeffcients of variation
        cv_c(:) = 0d0
        cv_n(:) = 0d0
        cv_h(:) = 0d0
        cv_l(:) = 0d0
        corr_hn(:) = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        cv_c(ij) = cv_c(ij) + c(ij, ia, ip, is)**2*phi(ij, ia, ip, is)
                        cv_n(ij) = cv_n(ij) + l(ij, ia, ip, is)**2*phi(ij, ia, ip, is)
                        cv_l(ij) = cv_l(ij) + (eff(ij)*theta(ip)*eta(is)* &
                                               l(ij, ia, ip, is))**2*phi(ij, ia, ip, is)
                        cv_h(ij) = cv_h(ij) + (eff(ij)*theta(ip)*eta(is))**2*phi(ij, ia, ip, is)
                        corr_hn(ij) = corr_hn(ij) + eff(ij)*theta(ip)*eta(is)* &
                                      l(ij, ia, ip, is)*phi(ij, ia, ip, is)
                    enddo
                enddo
            enddo
        enddo

        corr_hn = (corr_hn - n_coh*h_coh)/max(sqrt(cv_n - n_coh**2)*sqrt(cv_h - h_coh**2), 1d-10)

        cv_c = sqrt(cv_c - c_coh**2)/c_coh
        cv_n = sqrt(cv_n - n_coh**2)/max(n_coh, 1d-10)
        cv_h = sqrt(cv_h - h_coh**2)/max(h_coh, 1d-10)
        cv_l = sqrt(cv_l - l_coh**2)/max(l_coh, 1d-10)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none

        integer :: ij, ia, ip, is, iamax(JJ), ages(JJ)
        real*8 :: frac_bor(JJ)

        ! check for the maximium grid point used
        call check_grid(iamax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        write(21,'(a,a)')' IJ      CONS     HOURS  EARNINGS    INCOME      PENS    ASSETS', &
            '     CV(C)     CV(N)     CV(L)     VALUE     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,i10)')ij, c_coh(ij), n_coh(ij), l_coh(ij), &
                w*l_coh(ij) + r*a_coh(ij), pen(ij), a_coh(ij), &
                cv_c(ij), cv_n(ij), cv_l(ij), v_coh(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption')
        call plot(dble(ages), l_coh+pen, legend='Earnings')
        call plot(dble(ages), n_coh, legend='Hours')
        call execplot(xlabel='Age j', ylabel='Mean')

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

        call plot(dble(ages), cv_c, legend='Consumption')
        call plot(dble(ages), cv_l, legend='Earnings')
        call plot(dble(ages), cv_n, legend='Hours')
        call execplot(xlabel='Age j', ylabel='Coefficient of Variation')

        call plot(dble(ages), cv_l, legend='Earnings')
        call plot(dble(ages), cv_n, legend='Hours')
        call plot(dble(ages), cv_h, legend='Productivity')
        call plot(dble(ages), corr_hn, legend='Correlation')
        call execplot(xlabel='Age j', ylabel='Variance Decomposition Earnings')

        ! fraction of borrowing constrained households
        frac_bor = 0d0
        do ij = 1, JJ - 1
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        if(aplus(ij, ia, ip, is) < 1d-6)then
                            frac_bor(ij) = frac_bor(ij) + phi(ij, ia, ip, is)
                        endif
                    enddo
                enddo
            enddo
        enddo
        frac_bor(JJ) = 1d0

        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1))
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained Households')

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
