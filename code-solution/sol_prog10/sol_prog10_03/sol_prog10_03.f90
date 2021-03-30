!##############################################################################
! PROGRAM HyperbolicDiscounting
!
! ## The baseline life cycle model with hyperbolic discounting
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog10_03m.f90"

program HyperbolicDiscounting

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

        ! old-age transfers
        pen = 0d0
        pen(JR:JJ) = 0.5d0*w*sum(eff)/dble(JR-1)

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
    subroutine solve_household()

        implicit none
        integer :: ij, ia, ip, ip_max, is, is_max, ih
        real*8 :: x_in, wage
        logical :: check

        ! get decision in the last period of life
        do ih = 0, 1
            do ia = 0, NA
                aplus(JJ, ia, :, :, ih) = 0d0
                c(JJ, ia, :, :, ih) = (1d0+r)*a(ia) + pen(JJ)
                V(JJ, ia, :, :, ih) = valuefunc(0d0, c(JJ, ia, 1, 1, ih), JJ, 1, 1, ih)
            enddo
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

            do ih = 0, 1
                do ia = 0, NA

                    ! determine decision for zero assets at retirement without pension
                    if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                        aplus(ij, ia, :, :, ih) = 0d0
                        c(ij, ia, :, :, ih) = 0d0
                        V(ij, ia, :, :, ih) = valuefunc(0d0, 0d0, ij, 1, 1, ih)
                        cycle
                    endif

                    do ip = 1, ip_max
                        do is = 1, is_max

                            ! get initial guess for the individual choices
                            x_in = aplus(ij+1, ia, ip, is, ih)

                            ! set up communication variables
                            ij_com = ij
                            ia_com = ia
                            ip_com = ip
                            is_com = is

                            ! set hyperbolic discount factor
                            if(ih == 0)then
                                hdelta_com = 1d0
                            else
                                hdelta_com = hdelta
                            endif

                            ! solve the household problem using rootfinding
                            call fzero(x_in, foc, check)

                            ! write screen output in case of a problem
                            if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, is, ih

                            ! check for borrowing constraint
                            if(x_in < 0d0)then
                                x_in = 0d0
                                wage = eff(ij)*theta(ip)*eta(is)
                                cons_com = (1d0 + r)*a(ia) + wage + pen(ij)
                            endif

                            ! copy decisions
                            aplus(ij, ia, ip, is, ih) = x_in
                            c(ij, ia, ip, is, ih) = cons_com
                            V(ij, ia, ip, is, ih) = valuefunc(x_in, cons_com, ij, ip, is, ih)
                        enddo

                        ! copy decision in retirement age
                        if(ij >= JR)then
                            aplus(ij, ia, :, :, ih) = aplus(ij, ia, 1, 1, ih)
                            c(ij, ia, :, :, ih) = c(ij, ia, 1, 1, ih)
                            V(ij, ia, :, :, ih) = V(ij, ia, 1, 1, ih)
                        endif
                    enddo
                enddo
            enddo

            ! interpolate individual RHS
            call interpolate(ij)

            write(*,'(a,i3,a)')'Age: ', ij, ' DONE!'
        enddo

    end subroutine

    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, ip, is, is_p, ih
        real*8 :: chelp

        ! get RHS of first order condition
        do ia = 0, NA
            do ip = 1, NP
                do is = 1, NS
                    RHS(ij, ia, ip, is) = 0d0
                    do is_p = 1, NS
                        chelp = max(c(ij, ia, ip, is_p, 0), 1d-10)
                        RHS(ij, ia, ip, is) = RHS(ij, ia, ip, is) + &
                                              pi(is, is_p)*margu(chelp)
                    enddo
                    RHS(ij, ia, ip, is) = (beta*psi(ij)*(1d0+r)* &
                                                    RHS(ij, ia, ip, is))**(-gamma)
                enddo
            enddo
        enddo

        ! get expected value functions
        do ih = 0, 1
            do ia = 0, NA
                do ip = 1, NP
                    do is = 1, NS
                        EV(ij, ia, ip, is, ih) = 0d0
                        do is_p = 1, NS
                            EV(ij, ia, ip, is, ih) = EV(ij, ia, ip, is, ih) + &
                                                     pi(is, is_p)*V(ij, ia, ip, is_p, ih)
                        enddo
                        EV(ij, ia, ip, is, ih) = (egam*EV(ij, ia, ip, is, ih))**(1d0/egam)
                    enddo
                enddo
            enddo
        enddo

    end subroutine

    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, ip, is, is_p, ih, ial, iar
        real*8 :: varphi

        do ih = 0, 1

            ! set distribution to zero
            phi(:, :, :, :, ih) = 0d0

            ! get initial distribution in age 1
            do ip = 1, NP
                phi(1, 0, ip, is_initial, ih) = dist_theta(ip)
            enddo

            ! successively compute distribution over ages
            do ij = 2, JJ

                ! iterate over yesterdays gridpoints
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS

                            ! interpolate yesterday's savings decision
                            call linint_Grow(aplus(ij-1, ia, ip, is, ih), a_l, a_u, a_grow, NA, ial, iar, varphi)

                            ! restrict values to grid just in case
                            ial = min(ial, NA)
                            iar = min(iar, NA)
                            varphi = min(varphi, 1d0)

                            ! redistribute households
                            do is_p = 1, NS
                                phi(ij, ial, ip, is_p, ih) = phi(ij, ial, ip, is_p, ih) + &
                                                             pi(is, is_p)*varphi*phi(ij-1, ia, ip, is, ih)
                                phi(ij, iar, ip, is_p, ih) = phi(ij, iar, ip, is_p, ih) + &
                                                             pi(is, is_p)*(1d0-varphi)*phi(ij-1, ia, ip, is, ih)
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
        integer :: ij, ia, ip, is, ih

        do ih = 0, 1

            ! calculate cohort averages
            c_coh(:, ih) = 0d0
            y_coh(:, ih) = 0d0
            a_coh(:, ih) = 0d0
            v_coh(:, ih) = 0d0
            do ij = 1, JJ
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS
                            c_coh(ij, ih) = c_coh(ij, ih) + c(ij, ia, ip, is, ih)*phi(ij, ia, ip, is, ih)
                            y_coh(ij, ih) = y_coh(ij, ih) + w*eff(ij)*theta(ip)*eta(is)*phi(ij, ia, ip, is, ih)
                            a_coh(ij, ih) = a_coh(ij, ih) + a(ia)*phi(ij, ia, ip, is, ih)
                            v_coh(ij, ih) = v_coh(ij, ih) + V(ij, ia, ip, is, ih)*phi(ij, ia, ip, is, ih)
                        enddo
                    enddo
                enddo
            enddo

            ! calculate cohort specific coeffcients of variation
            cv_c(:, ih) = 0d0
            cv_y(:, ih) = 0d0
            do ij = 1, JJ
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS
                            cv_c(ij, ih) = cv_c(ij, ih) + c(ij, ia, ip, is, ih)**2*phi(ij, ia, ip, is, ih)
                            cv_y(ij, ih) = cv_y(ij, ih) + (w*eff(ij)*theta(ip)*eta(is))**2*phi(ij, ia, ip, is, ih)
                        enddo
                    enddo
                enddo
            enddo
        enddo
        cv_c = sqrt(cv_c - c_coh**2)/c_coh
        cv_y = sqrt(cv_y - y_coh**2)/max(y_coh, 1d-10)

    end subroutine

    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ih, ij, ia, ip, is, iamax(JJ, 0:1), ages(JJ)
        real*8 :: frac_bor(JJ, 0:1)

        ! check for the maximium grid point used
        call check_grid(iamax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        ! set up output file
        do ih = 0, 1

            if(ih == 0)then
                write(21,'(a/)')'RATIONAL AGENT'
            else
                write(21,'(a/)')'HYPERBOLIC DISCOUNTING AGENT'
            endif

            write(21,'(a,a)')' IJ      CONS  EARNINGS    INCOME      PENS    ASSETS', &
                '     CV(C)     CV(L)     VALUE     IAMAX'
            do ij = 1, JJ
                write(21,'(i3,8f10.3,i10)')ij, c_coh(ij, ih), y_coh(ij, ih), &
                    y_coh(ij, ih)+r*a_coh(ij, ih), pen(ij), a_coh(ij, ih), &
                    cv_c(ij, ih), cv_y(ij, ih), v_coh(ij, ih), iamax(ij, ih)
            enddo
            write(21,'(a/)')'--------------------------------------------------------------------'

        enddo

        ! plot output
        call plot(dble(ages), c_coh(:, 0), legend='Rational')
        call plot(dble(ages), c_coh(:, 1), legend='Hyperbolic')
        call plot(dble(ages), y_coh(:, 0)+pen, legend='Earnings')
        call plot(dble(ages), y_coh(:, 1)+pen, legend='Earnings')
        call execplot(xlabel='Age j', ylabel='Consumption')

        call plot(dble(ages), a_coh(:, 0), legend='Rational')
        call plot(dble(ages), a_coh(:, 1), legend='Hyperbolic')
        call execplot(xlabel='Age j', ylabel='Assets')

        call plot(dble(ages), cv_c(:, 0), legend='Rational')
        call plot(dble(ages), cv_c(:, 1), legend='Hyperbolic')
        call execplot(xlabel='Age j', ylabel='Coefficient of Variation')

        ! calculate fraction of borrowing constrained households
        frac_bor = 0d0
        do ih = 0, 1
            do ij = 1, JJ-1
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS
                            if(aplus(ij, ia, ip, is, ih) < 1d-6)then
                                frac_bor(ij, ih) = frac_bor(ij, ih) + phi(ij, ia, ip, is, ih)
                            endif
                        enddo
                    enddo
                enddo
            enddo
            frac_bor(JJ, ih) = 1d0
        enddo

        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 0), legend='Rational')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 1), legend='Hyperbolic')
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained Households')

    end subroutine

    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax)

        implicit none
        integer :: iamax(JJ, 0:1), ij, ia, ip, is, ih

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ih = 0, 1
                do ia = 0, NA
                    do ip = 1, NP
                        do is = 1, NS
                            if(phi(ij, ia, ip, is, ih) > 1d-8)iamax(ij, ih) = ia
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
