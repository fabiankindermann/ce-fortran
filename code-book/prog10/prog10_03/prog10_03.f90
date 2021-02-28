!##############################################################################
! PROGRAM FemaleParticipation
!
! ## A life cycle model with children and female labor force participation
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
include "prog10_03m.f90"

program FemaleParticipation

    use globals

    implicit none

    ! initialize remaining variables
    call initialize()

    ! start the clock
    call tic()

    ! solve the household problem
    call solve_household()

    ! calculate the distribution of households over state space
    call get_distribution()

    ! aggregate individual decisions
    call aggregation()

    ! stop the clock
    call toc()

    call output()

    ! close files
    close(21)

contains


    ! initializes all remaining variables
    subroutine initialize

        implicit none
        real*8 :: price
        integer :: children(10), ic, ij

        ! net prices (after taxes and tranfers)
        r = 0.04d0
        w_m = 1.00d0
        w_f = 0.75d0

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

        ! old-age transfers
        pen = 0d0
        pen(JR:JJ) = 0.8d0*sum(eff)/dble(JR-1)

        ! initialize fixed effect
        dist_theta = 1d0/dble(NP)
        theta(1)   = exp(-sqrt(sigma_theta))
        theta(2)   = exp( sqrt(sigma_theta))

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize the human capital grid by defining human capital maximum
        h_l = 0d0
        h_max(1) = 0d0
        do ij = 2, JR-1
            h_max(ij) = h_max(ij-1) + xi(1) + xi(2)*(ij-1)
        enddo
        h_max(JR:JJ) = h_max(JR-1)
        h_u = maxval(h_max)

        call grid_Cons_Grow(h, h_l, h_u, h_grow)

        ! define arrival of children
        children = 0
        children(1:2) = (/30, 32/)


        ! different setups for Figure 10.5 (uncomment respective line)
        !children = 0                     ! No Children
        !children(1:2) = (/23, 25/)       ! Children 23/25
        !children(1:3) = (/27, 30, 33/)   ! Children 27/30/33

        ! set respective child variables
        nchild = 0
        pchild = 0d0
        price = w_f*exp(h_max(10))
        do ic = 1, 10
            if(children(ic) > 20 .and. children(ic) < 50)then

                ! calculate age child is born
                ij = children(ic) - 20

                ! set number of children
                nchild(ij:ij+17) = nchild(ij:ij+17) + 1

                ! set cost of child care
                pchild(ij:ij+2)     = pchild(ij:ij+2)    + 1.00d0*price
                pchild(ij+3:ij+5)   = pchild(ij+3:ij+5)  + 0.80d0*price
                pchild(ij+6:ij+11)  = pchild(ij+6:ij+11) + 0.60d0*price
                pchild(ij+12:ij+17) = pchild(ij+12:ij+17)+ 0.40d0*price
            endif
        enddo

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, ih, ih_max, ip, ip_max, ism, isf, is_max, il, il_max
        real*8 :: x_in, wage_m, wage_f
        real*8 :: c_temp(0:1), ap_temp(0:1), u_temp(0:1)
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            aplus(JJ, ia, :, :, :, :) = 0d0
            c(JJ, ia, :, :, :, :) = (1d0+r)*a(ia) + pen(JJ)
            l(JJ, ia, :, :, :, :) = 0d0
            V(JJ, ia, :, :, :, :) = valuefunc(0d0, c(JJ, ia, 0, 1, 1, 1), 0, JJ, 0, 1, 1, 1)
        enddo

        ! interpolate individual RHS
        call interpolate(JJ)

        do ij = JJ-1, 1, -1

            ! check about how many states to iterate
            if(ij >= JR)then
                ih_max = 0
                ip_max = 1
                is_max = 1
                il_max = 0
                u_temp(1) = -1d100
            else
                ih_max = NH
                ip_max = NP
                is_max = NS
                il_max = 1
            endif

            do ia = 0, NA

                ! determine decision for zero assets at retirement without pension
                if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                    aplus(ij, ia, :, :, :, :) = 0d0
                    c(ij, ia, :, :, :, :) = 0d0
                    l(ij, ia, :, :, :, :) = 0d0
                    V(ij, ia, :, :, :, :) = valuefunc(0d0, 0d0, 0, ij, 0, 1, 1, 1)
                    cycle
                endif

                do ih = 0, ih_max

                    ! check whether h(ih) is greater than h_max(ij)
                    if(h(ih) > h_max(ij))then
                        aplus(ij, ia, ih, :, :, :) = aplus(ij, ia, ih-1, :, :, :)
                        c(ij, ia, ih, :, :, :) = c(ij, ia, ih-1, :, :, :)
                        l(ij, ia, ih, :, :, :) = l(ij, ia, ih-1, :, :, :)
                        V(ij, ia, ih, :, :, :) = V(ij, ia, ih-1, :, :, :)
                        cycle
                    endif

                    do ip = 1, ip_max
                        do ism = 1, is_max
                            do isf = 1, is_max

                                ! determine solution for both working decisions
                                do il = 0, il_max

                                    ! get initial guess for the individual choices
                                    x_in = aplus(ij+1, ia, ih, ip, ism, isf)

                                    ! set up communication variables
                                    ij_com = ij
                                    ia_com = ia
                                    ih_com = ih
                                    ip_com = ip
                                    ism_com = ism
                                    isf_com = isf
                                    il_com = il

                                    ! solve the household problem using rootfinding
                                    call fzero(x_in, foc, check)

                                    ! check for borrowing constraint
                                    if(x_in < 0d0)then
                                        x_in = 0d0
                                        wage_m = w_m*eff(ij)*theta(ip)*eta(ism)
                                        wage_f = w_f*dble(il)*(exp(h(ih))*theta(ip)*eta(isf) - pchild(ij))
                                        cons_com = (1d0+r)*a(ia) + pen(ij) + wage_f + wage_m
                                    endif

                                    ! write screen output in case of a problem
                                    if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ia, ih, ip, ism, isf

                                    c_temp(il) = cons_com
                                    ap_temp(il) = x_in
                                    u_temp(il) = valuefunc(x_in, cons_com, il, ij, ih, ip, ism, isf)
                                enddo

                                ! choose the labor force status the gives more utility
                                if(u_temp(1) >= u_temp(0))then
                                    aplus(ij, ia, ih, ip, ism, isf) = ap_temp(1)
                                    c(ij, ia, ih, ip, ism, isf) = c_temp(1)
                                    l(ij, ia, ih, ip, ism, isf) = 1d0
                                    V(ij, ia, ih, ip, ism, isf) = u_temp(1)
                                else
                                    aplus(ij, ia, ih, ip, ism, isf) = ap_temp(0)
                                    c(ij, ia, ih, ip, ism, isf) = c_temp(0)
                                    l(ij, ia, ih, ip, ism, isf) = 0d0
                                    V(ij, ia, ih, ip, ism, isf) = u_temp(0)
                                endif
                            enddo
                        enddo
                    enddo
                enddo

                ! copy decision in retirement age
                if(ij >= JR)then
                    aplus(ij, ia, :, :, :, :) = aplus(ij, ia, 0, 1, 1, 1)
                    c(ij, ia, :, :, :, :) = c(ij, ia, 0, 1, 1, 1)
                    l(ij, ia, :, :, :, :) = l(ij, ia, 0, 1, 1, 1)
                    V(ij, ia, :, :, :, :) = V(ij, ia, 0, 1, 1, 1)
                endif
            enddo

            ! interpolate individual RHS
            call interpolate(ij)

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, ih, ip, ism, isf, ism_p, isf_p, ihl, ihr, il
        integer :: ih_max, ip_max, is_max, il_max
        real*8 :: chelp, Vhelp, hplus, varphi

        do ia = 0, NA

            ! check about how many is to iterate
            if(ij > JR)then
                ih_max = 0
                ip_max = 1
                is_max = 1
                il_max = 0
            else
                ih_max = NH
                ip_max = NP
                is_max = NS
                il_max = 1
            endif

            do ih = 0, ih_max
                do ip = 1, ip_max
                    do ism = 1, is_max
                        do isf = 1, is_max
                            do il = 0, il_max
                                RHS(ia, ih, ip, ism, isf, il) = 0d0
                                EV(ia, ih, ip, ism, isf, il) = 0d0

                                ! interpolate human capital for tomorrow
                                hplus = h(ih) + (xi(1) + xi(2)*dble(ij-1))*dble(il) - del_h*dble(1-il)
                                hplus = max(hplus, h_l)
                                call linint_Grow(hplus, h_l, h_u, h_grow, NH, ihl, ihr, varphi)

                                ! iterate over all potential future states
                                do ism_p = 1, NS
                                    do isf_p = 1, NS

                                        ! right hand side of first order condition
                                        chelp = varphi      *c(ij, ia, ihl, ip, ism_p, isf_p) + &
                                                (1d0-varphi)*c(ij, ia, ihr, ip, ism_p, isf_p)
                                        chelp = max(chelp, 1d-10)
                                        RHS(ia, ih, ip, ism, isf, il) = RHS(ia, ih, ip, ism, isf, il) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)*margu(chelp, ij)

                                        ! expected value function
                                        Vhelp = varphi      *(egam*V(ij, ia, ihl, ip, ism_p, isf_p))**(1d0/egam)+ &
                                                (1d0-varphi)*(egam*V(ij, ia, ihr, ip, ism_p, isf_p))**(1d0/egam)
                                        Vhelp = max(Vhelp, 1d-10)**egam/egam
                                        EV(ia, ih, ip, ism, isf, il) = EV(ia, ih, ip, ism, isf, il) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)*Vhelp
                                    enddo
                                enddo
                                RHS(ia, ih, ip, ism, isf, il) = &
                                        ((1d0+r)*beta*psi(ij)*RHS(ia, ih, ip, ism, isf, il))**(-gamma)
                                EV(ia, ih, ip, ism, isf, il) = &
                                        (egam*EV(ia, ih, ip, ism, isf, il))**(1d0/egam)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, ih, ip, ism, ism_p, isf, isf_p
        integer :: ial, iar, ihl, ihr
        real*8 :: varphi_a, varphi_h, htoday, labor

        ! set distribution to zero
        phi(:, :, :, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ip = 1, NP
            phi(1, 0, 0, ip, is_initial, is_initial) = dist_theta(ip)
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do ih = 0, NH
                    do ip = 1 , NP
                        do ism = 1, NS
                            do isf = 1, NS

                                ! interpolate yesterday's savings decision
                                call linint_Grow(aplus(ij-1, ia, ih, ip, ism, isf), a_l, a_u, &
                                        a_grow, NA, ial, iar, varphi_a)

                                ! restrict values to grid just in case
                                ial = min(ial, NA)
                                iar = min(iar, NA)
                                varphi_a = min(varphi_a, 1d0)

                                ! get today's human capital
                                labor = l(ij-1, ia, ih, ip, ism, isf)
                                htoday = h(ih) + (xi(1) + xi(2)*dble(ij-1))*labor - del_h*(1d0-labor)
                                htoday = max(htoday, h_l)
                                call linint_Grow(htoday, h_l, h_u, h_grow, NH, ihl, ihr, varphi_h)
                                ihl = min(ihl, NA)
                                ihr = min(ihr, NA)
                                varphi_h = min(varphi_h, 1d0)

                                if(ij >= JR)then
                                    ihl = ih
                                    ihr = ih
                                endif

                                ! redistribute households
                                do ism_p = 1, NS
                                    do isf_p = 1, NS
                                        phi(ij, ial, ihl, ip, ism_p, isf_p) = &
                                                phi(ij, ial, ihl, ip, ism_p, isf_p) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)* &
                                                varphi_a*varphi_h*phi(ij-1, ia, ih, ip, ism, isf)
                                        phi(ij, iar, ihl, ip, ism_p, isf_p) = &
                                                phi(ij, iar, ihl, ip, ism_p, isf_p) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)* &
                                                (1d0-varphi_a)*varphi_h*phi(ij-1, ia, ih, ip, ism, isf)
                                        phi(ij, ial, ihr, ip, ism_p, isf_p) = &
                                                phi(ij, ial, ihr, ip, ism_p, isf_p) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)* &
                                                varphi_a*(1d0-varphi_h)*phi(ij-1, ia, ih, ip, ism, isf)
                                        phi(ij, iar, ihr, ip, ism_p, isf_p) = &
                                                phi(ij, iar, ihr, ip, ism_p, isf_p) + &
                                                pi(ism, ism_p)*pi(isf, isf_p)* &
                                                (1d0-varphi_a)*(1d0-varphi_h)*phi(ij-1, ia, ih, ip, ism, isf)
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ia, ih, ip, ism, isf
        real*8 :: wage_m, wage_f

        ! calculate cohort averages
        c_coh(:)  = 0d0
        ym_coh(:) = 0d0
        yf_coh(:) = 0d0
        l_coh(:)  = 0d0
        a_coh(:)  = 0d0
        h_coh(:)  = 0d0
        v_coh(:)  = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do ih = 0, NH
                    do ip = 1, NP
                        do ism = 1, NS
                            do isf = 1, NS
                                wage_m = w_m*eff(ij)*theta(ip)*eta(ism)
                                wage_f = w_f*exp(h(ih))*theta(ip)*eta(isf)*l(ij, ia, ih, ip, ism, isf)
                                c_coh(ij)  = c_coh(ij) + c(ij, ia, ih, ip, ism, isf)*phi(ij, ia, ih, ip, ism, isf)
                                ym_coh(ij) = ym_coh(ij)+ wage_m*phi(ij, ia, ih, ip, ism, isf)
                                yf_coh(ij) = yf_coh(ij)+ wage_f*phi(ij, ia, ih, ip, ism, isf)
                                l_coh(ij)  = l_coh(ij) + l(ij, ia, ih, ip, ism, isf)*phi(ij, ia, ih, ip, ism, isf)
                                a_coh(ij)  = a_coh(ij) + a(ia)*phi(ij, ia, ih, ip, ism, isf)
                                h_coh(ij)  = h_coh(ij) + exp(h(ih))*phi(ij, ia, ih, ip, ism, isf)
                                v_coh(ij)  = v_coh(ij) + V(ij, ia, ih, ip, ism, isf)*phi(ij, ia, ih, ip, ism, isf)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! calculate cohort specific variances of logs
        cv_c = 0d0
        cv_y = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do ih = 0, NH
                    do ip = 1, NP
                        do ism = 1, NS
                            do isf = 1, NS
                                wage_m = w_m*eff(ij)*theta(ip)*eta(ism)
                                wage_f = w_f*exp(h(ih))*theta(ip)*eta(isf)*l(ij, ia, ih, ip, ism, isf)
                                cv_c(ij) = cv_c(ij) + c(ij, ia, ih, ip, ism, isf)**2*phi(ij, ia, ih, ip, ism, isf)
                                cv_y(ij) = cv_y(ij) + (wage_m+wage_f)**2*phi(ij, ia, ih, ip, ism, isf)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        cv_c = sqrt(cv_c-c_coh**2)/c_coh
        cv_y = sqrt(cv_y-(ym_coh+yf_coh)**2)/max(ym_coh+yf_coh, 1d-10)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ia, ih, ip, ism, isf, iamax(JJ), ages(JJ)
        real*8 :: frac_bor(JJ)

        ! check for the maximium grid point used
        call check_grid(iamax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        write(21, '(a,a)')' IJ      CONS      FLFP    HUMCAP    WAGE_F    WAGE_M  EARNINGS', &
            '      PENS    ASSETS     CV(C)     CV(Y)     VALUE     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,11f10.3,i10)')ij, c_coh(ij), l_coh(ij), h_coh(ij), &
                    yf_coh(ij), ym_coh(ij), yf_coh(ij)+ym_coh(ij), pen(ij), a_coh(ij), &
                    cv_c(ij), cv_y(ij), v_coh(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption')
        call plot(dble(ages), c_coh/sqrt(2d0+nchild), legend='Normalized Cons.')
        call plot(dble(ages), ym_coh+yf_coh+pen, legend='Total Earnings')
        call execplot(xlabel='Age j', ylabel='Mean', ylim = (/0d0, 5d0/))

        call plot(dble(ages), ym_coh, legend='Male')
        call plot(dble(ages), yf_coh, legend='Female')
        call plot(dble(ages), ym_coh+yf_coh, legend='Total')
        call execplot(xlabel='Age j', ylabel='Earnings', ylim = (/0d0, 5d0/))

        call plot(dble(ages), l_coh)
        call execplot(xlabel='Age j', ylabel='Female Labor Force Participation', ylim = (/0d0, 1d0/))

        call plot(dble(ages(1:JR-1)), eff(1:JR-1), legend='Men')
        call plot(dble(ages(1:JR-1)), h_coh(1:JR-1), legend='Women')
        call execplot(xlabel='Age j', ylabel='Human Capital', ylim = (/0d0, 2d0/))

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

        ! fraction of borrowing constrained households
        frac_bor = 0d0
        do ij = 1, JJ-1
            do ia = 0, NA
                do ih = 0, NH
                    do ip = 1, NP
                        do ism = 1, NS
                            do isf = 1, NS
                                if(aplus(ij, ia, ih, ip, ism, isf) < 1d-6)then
                                    frac_bor(ij) = frac_bor(ij) + phi(ij, ia, ih, ip, ism, isf)
                                endif
                            enddo
                        enddo
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
        integer :: iamax(JJ), ij, ia, ih, ip, ism, isf

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do ih = 0, NH
                    do ip = 1, NP
                        do ism = 1, NS
                            do isf = 1, NS
                                if(phi(ij, ia, ih, ip, ism, isf) > 1d-8)iamax(ij) = ia
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
