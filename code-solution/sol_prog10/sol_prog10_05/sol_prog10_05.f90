!##############################################################################
! PROGRAM HealthRiskCost
!
! ## The baseline life cycle model with health risk
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
include "sol_prog10_05m.f90"

program HealthRiskCost

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
        integer :: im, ip, ij, ijj
        real*8 :: abor_temp, k(JJ)

        ! net prices (after taxes and tranfers)
        r = 0.04d0
        w = 1.0d0

        ! set survival probabilities
        psi(:, 0) = (/1.00000d0, 0.99923d0, 0.99914d0, 0.99914d0, 0.99912d0, &
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

        ! health depended survival probabilities
        psi(:, 1) = chi*psi(:, 0)

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

        ! initialize health shock
        dist_m(0) = 0.95d0
        dist_m(1) = 0.05d0

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! probability to have bad health when current health is good
        call grid_Cons_Equi(pi_m(:, 1, 0, 1), 0.1d0, 0.4d0)
        call grid_Cons_Equi(pi_m(:, 2, 0, 1), 0.0d0, 0.2d0)
        pi_m(:, :, 0, 0) = 1d0 - pi_m(:, :, 0, 1)

        ! probability to have bad health when current health is bad
        call grid_Cons_Equi(pi_m(:, 1, 1, 1), 0.6d0, 0.9d0)
        call grid_Cons_Equi(pi_m(:, 2, 1, 1), 0.6d0, 0.9d0)
        pi_m(:, :, 1, 0) = 1d0 - pi_m(:, :, 1, 1)

        ! initialize impact of shock on productivity
        varrho(1) = 0.2d0
        varrho(2) = 0.1d0
        do ip = 1, NP
            do im = 0, NM
                varrho_m(ip, im) = exp(-varrho(ip)*dble(im))
            enddo
        enddo

        ! calculate the out-of-pocket expenses for health
        call grid_Cons_Equi(k, 0.3d0, 0.9d0)

        ! normally distributed stochastic term
        call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
        zeta = exp(zeta)

        ! out of pocket expenses
        do ij = 1, JJ
            hc(ij, :) = k(ij)*zeta
        enddo

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! value of the borrowing constraint
        a_bor = 0d0

        ! calculate endogenous borrowing constraints
        do ij = 1, JJ
            do ip = 1, NP

                ! calculate natural borrowing limit
                abor_temp = 0d0
                do ijj = JJ, ij, -1
                    abor_temp = abor_temp/(1d0+r) + eff(ijj)*theta(ip)*eta(1) + pen(ijj)
                enddo
                abor_temp = min(-abor_temp/(1d0+r)+1d-4, 0d0)

                ! set maximum of natural and exogenous borrowing constraint
                a_bor(ij, ip) = max(a_bor(ij, ip), abor_temp)
            enddo
        enddo

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, ip, ip_max, im, is, is_max, ih
        real*8 :: x_in, wage
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            aplus(JJ, ia, :, :, :, :) = 0d0
            do ih = 1, NH
                c(JJ, ia, :, :, :, ih) = (1d0+r)*a(ia) + pen(JJ) - hc(JJ, ih) + &
                                         max(c_floor + hc(JJ, ih) - (1d0+r)*a(ia) - pen(JJ), 0d0)
                do im = 0, NM
                    V(JJ, ia, :, im, :, ih) = valuefunc(0d0, c(JJ, ia, 1, im, 1, ih), JJ, 1, im, 1)
                enddo
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

            do ia = 0, NA

                ! determine decision for zero assets at retirement without pension
                if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                    aplus(ij, ia, :, :, :, :) = 0d0
                    c(ij, ia, :, :, :, :) = 0d0
                    do im = 0, NM
                        V(ij, ia, :, im, :, :) = valuefunc(0d0, 0d0, ij, 1, im, 1)
                    enddo
                    cycle
                endif

                ! solve household problem using rootfinding method
                do ip = 1, ip_max
                    do im = 0, NM
                        do is = 1, is_max
                            do ih = 1, NH

                                ! get initial guess for the individual choices
                                x_in = aplus(ij+1, ia, ip, im, is, ih)

                                ! set up communication variables
                                ij_com = ij
                                ia_com = ia
                                ip_com = ip
                                im_com = im
                                is_com = is
                                ih_com = ih

                                ! solve the household problem using rootfinding
                                call fzero(x_in, foc, check)

                                ! write screen output in case of a problem
                                if(check)write(*,'(a, 6i4)')'ERROR IN ROOTFINDING : ', ij, ia, ip, im, is, ih

                                ! check for borrowing constraint
                                if(x_in < 0d0)then
                                    x_in = 0d0
                                    wage = w*eff(ij)*theta(ip)*eta(is)*varrho_m(ip, im)
                                    cons_com = max((1d0+r)*a(ia) + wage + pen(ij) - hc(ij, ih) &
                                    + max(c_floor + hc(ij, ih) - (1d0+r)*a(ia) - wage - pen(ij), 0d0), 1d-10)
                                endif

                                ! copy decisions
                                aplus(ij, ia, ip, im, is, ih) = x_in
                                c(ij, ia, ip, im, is, ih) = cons_com
                                V(ij, ia, ip, im, is, ih) = valuefunc(x_in, cons_com, ij, ip, im, is)
                            enddo
                        enddo
                    enddo
                enddo

                ! copy decision in retirement age
                if(ij >= JR)then
                    do im = 0, NM
                        do ih = 1, NH
                            aplus(ij, ia, :, im, :, ih) = aplus(ij, ia, 1, im, 1, ih)
                            c(ij, ia, :, im, :, ih) = c(ij, ia, 1, im, 1, ih)
                            V(ij, ia, :, im, :, ih) = V(ij, ia, 1, im, 1, ih)
                        enddo
                    enddo
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
        integer :: ia, ip, im, im_p, is, is_p, ih_p
        real*8 :: chelp

        do ia = 0, NA
            do ip = 1, NP
                do im = 0, NM
                    do is = 1, NS

                        ! calculate the RHS of the first order condition
                        RHS(ij, ia, ip, im, is) = 0d0
                        EV(ij, ia, ip, im, is) = 0d0
                        do im_p = 0, NM
                            do is_p = 1, NS
                                do ih_p = 1, NH
                                    chelp = max(c(ij, ia, ip, im_p, is_p, ih_p), 1d-10)
                                    RHS(ij, ia, ip, im, is) = RHS(ij, ia, ip, im, is) + dist_zeta(ih_p)* &
                                                pi_m(ij, ip, im, im_p)*pi(is, is_p)*margu(chelp, im_p)
                                    EV(ij, ia, ip, im, is) = EV(ij, ia, ip, im, is) + dist_zeta(ih_p)* &
                                                pi_m(ij, ip, im, im_p)*pi(is, is_p)*V(ij, ia, ip, im_p, is_p, ih_p)
                                enddo
                            enddo
                        enddo
                        RHS(ij, ia, ip, im, is) = (beta*psi(ij, im)*(1d0+r)/(1d0-delta*dble(im)) &
                                                   *RHS(ij, ia, ip, im, is))**(-gamma)
                        EV(ij, ia, ip, im, is) = (egam*EV(ij, ia, ip, im, is))**(1d0/egam)
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, ip, im, im_p, is, is_p, ih, ih_p, ial, iar
        real*8 :: varphi

        ! set distribution to zero
        phi(:, :, :, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ip = 1, NP
            do im = 0, NM
                do ih = 1, NH

                    ! find the zero on the asset grid
                    call linint_Grow(-a_bor(1, ip), a_l, a_u, a_grow, NA, ial, iar, varphi)
                    phi(1, ial, ip, im, is_initial, ih) = varphi*dist_theta(ip)*dist_m(im)*dist_zeta(ih)
                    phi(1, iar, ip, im, is_initial, ih) = (1d0-varphi)*dist_theta(ip)*dist_m(im)*dist_zeta(ih)

                enddo
            enddo
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do ip = 1, NP
                    do im = 0, NM
                        do is = 1, NS
                            do ih = 1, NH

                                ! interpolate yesterday's savings decision
                                call linint_Grow(aplus(ij-1, ia, ip, im, is, ih), a_l, &
                                        a_u, a_grow, NA, ial, iar, varphi)

                                ! restrict values to grid just in case
                                ial = min(ial, NA)
                                iar = min(iar, NA)
                                varphi = min(varphi, 1d0)

                                ! redistribute households
                                do im_p = 0, NM
                                    do is_p = 1, NS
                                        do ih_p = 1, NH
                                            phi(ij, ial, ip, im_p, is_p, ih_p) = phi(ij, ial, ip, im_p, is_p, ih_p) + &
                                                 pi(is, is_p)*pi_m(ij, ip, im, im_p)*varphi*dist_zeta(ih_p)* &
                                                                                 phi(ij-1, ia, ip, im, is, ih)
                                            phi(ij, iar, ip, im_p, is_p, ih_p) = phi(ij, iar, ip, im_p, is_p, ih_p) + &
                                                 pi(is, is_p)*pi_m(ij, ip, im, im_p)*(1d0-varphi)*dist_zeta(ih_p)* &
                                                                                 phi(ij-1, ia, ip, im, is, ih)
                                        enddo
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
        integer :: ij, ia, ip, im, is, ih

        ! calculate fraction of good vs. bad health households
        do ij = 1, JJ
            do ip = 1, NP
                do im = 0, NM
                    frac_phi(ij, ip, im) = sum(phi(ij, :, ip, im, :, :))
                enddo
            enddo
        enddo

        ! calculate cohort averages
        c_coh(:, :, :) = 0d0
        y_coh(:, :, :) = 0d0
        a_coh(:, :, :) = 0d0
        v_coh(:, :, :) = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do ip = 1, NP
                    do im = 0, NM
                        do is = 1, NS
                            do ih = 1, NH
                                c_coh(ij, ip, im) = c_coh(ij, ip, im) + c(ij, ia, ip, im, is, ih) &
                                                    *phi(ij, ia, ip, im, is, ih)/frac_phi(ij, ip, im)
                                y_coh(ij, ip, im) = y_coh(ij, ip, im) + eff(ij)*theta(ip)*eta(is)**varrho_m(ip, im) &
                                                    *phi(ij, ia, ip, im, is, ih)/frac_phi(ij, ip, im)
                                a_coh(ij, ip, im) = a_coh(ij, ip, im) + a(ia) &
                                                    *phi(ij, ia, ip, im, is, ih)/frac_phi(ij, ip, im)
                                v_coh(ij, ip, im) = v_coh(ij, ip, im) + V(ij, ia, ip, im, is, ih) &
                                                    *phi(ij, ia, ip, im, is, ih)/frac_phi(ij, ip, im)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! recover unconditional cohort averages
        do ij = 1, JJ
            do ip = 1, NP
                do im = 0, NM
                    c_coh(ij, ip, NM+1) = c_coh(ij, ip, NM+1) + c_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, ip, :))
                    y_coh(ij, ip, NM+1) = y_coh(ij, ip, NM+1) + y_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, ip, :))
                    a_coh(ij, ip, NM+1) = a_coh(ij, ip, NM+1) + a_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, ip, :))
                    v_coh(ij, ip, NM+1) = v_coh(ij, ip, NM+1) + v_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, ip, :))
                    c_coh(ij, NP+1, im) = c_coh(ij, NP+1, im) + c_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, im))
                    y_coh(ij, NP+1, im) = y_coh(ij, NP+1, im) + y_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, im))
                    a_coh(ij, NP+1, im) = a_coh(ij, NP+1, im) + a_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, im))
                    v_coh(ij, NP+1, im) = v_coh(ij, NP+1, im) + v_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, im))
                    c_coh(ij, NP+1, NM+1) = c_coh(ij, NP+1, NM+1) + c_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, :))
                    y_coh(ij, NP+1, NM+1) = y_coh(ij, NP+1, NM+1) + y_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, :))
                    a_coh(ij, NP+1, NM+1) = a_coh(ij, NP+1, NM+1) + a_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, :))
                    v_coh(ij, NP+1, NM+1) = v_coh(ij, NP+1, NM+1) + v_coh(ij, ip, im) &
                                                *frac_phi(ij, ip, im)/sum(frac_phi(ij, :, :))
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ia, ip, im, is, ih, iamax(JJ), ages(JJ)
        real*8 :: frac_bor(JJ, NP+1, 0:NM+1)

        ! check for the maximium grid point used
        call check_grid(iamax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        ! set up output file
        do im = 0, NM + 1

            if(im == 0)then
                write(21,'(a/)')'GOOD HEALTH'
            elseif(im == 1)then
                write(21,'(a/)')'BAD HEALTH'
            else
                write(21,'(a/)')'UNCONDITIONAL COHORT AVERAGE'
            endif

            write(21,'(a,a)')' IJ      CONS  EARNINGS    INCOME      PENS    ASSETS', &
                '     VALUE     IAMAX'
            do ij = 1, JJ
                write(21,'(i3,6f10.3,i10)')ij, c_coh(ij, 3, im), y_coh(ij, 3, im), &
                    y_coh(ij, 3, im)+r*a_coh(ij, 3, im), pen(ij), a_coh(ij, 3, im), &
                    v_coh(ij, 3, im), iamax(ij)
            enddo
            write(21,'(a/)')'--------------------------------------------------------------------'

        enddo

        ! plot output
        call plot(dble(ages), c_coh(:, 1, 0), legend='Low Skilled - Good Health')
        call plot(dble(ages), c_coh(:, 1, 1), legend='Low Skilled - Bad Health')
        call plot(dble(ages), c_coh(:, 1, 2), legend='Low Skilled - Average')
        call plot(dble(ages), c_coh(:, 2, 0), legend='High Skilled - Good Health')
        call plot(dble(ages), c_coh(:, 2, 1), legend='High Skilled - Bad Health')
        call plot(dble(ages), c_coh(:, 2, 2), legend='High Skilled - Average')
        call plot(dble(ages), c_coh(:, 3, 0), legend='Good Health')
        call plot(dble(ages), c_coh(:, 3, 1), legend='Bad Health')
        call plot(dble(ages), c_coh(:, 3, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Consumption')

        call plot(dble(ages), y_coh(:, 1, 0)+pen, legend='Low Skilled - Good Health')
        call plot(dble(ages), y_coh(:, 1, 1)+pen, legend='Low Skilled - Bad Health')
        call plot(dble(ages), y_coh(:, 1, 2)+pen, legend='Low Skilled - Average')
        call plot(dble(ages), y_coh(:, 2, 0)+pen, legend='High Skilled - Good Health')
        call plot(dble(ages), y_coh(:, 2, 1)+pen, legend='High Skilled - Bad Health')
        call plot(dble(ages), y_coh(:, 2, 2)+pen, legend='High Skilled - Average')
        call plot(dble(ages), y_coh(:, 3, 0)+pen, legend='Good Health')
        call plot(dble(ages), y_coh(:, 3, 1)+pen, legend='Bad Health')
        call plot(dble(ages), y_coh(:, 3, 2)+pen, legend='Average')
        call execplot(xlabel='Age j', ylabel='Earnings')

        call plot(dble(ages), a_coh(:, 1, 0), legend='Low Skilled - Good Health')
        call plot(dble(ages), a_coh(:, 1, 1), legend='Low Skilled - Bad Health')
        call plot(dble(ages), a_coh(:, 1, 2), legend='Low Skilled - Average')
        call plot(dble(ages), a_coh(:, 2, 0), legend='High Skilled - Good Health')
        call plot(dble(ages), a_coh(:, 2, 1), legend='High Skilled - Bad Health')
        call plot(dble(ages), a_coh(:, 2, 2), legend='High Skilled - Average')
        call plot(dble(ages), a_coh(:, 3, 0), legend='Bad Health')
        call plot(dble(ages), a_coh(:, 3, 1), legend='Good Health')
        call plot(dble(ages), a_coh(:, 3, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Assets')

        ! calculate fraction of borrowing constrained households
        frac_bor = 0d0
        do ij = 1, JJ-1
            do ia = 0, NA
                do ip = 1, NP
                    do im = 0, NM
                        do is = 1, NS
                            do ih = 1, NH
                                if(aplus(ij, ia, ip, im, is, ih) < 1d-6)then
                                    frac_bor(ij, ip, im) = frac_bor(ij, ip, im) &
                                            + phi(ij, ia, ip, im, is, ih)/frac_phi(ij, ip, im)
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        frac_bor(JJ, :, :) = 1d0

        do ij = 1, JJ
            do ip = 1, NP
                do im = 0, NM
                    frac_bor(ij, ip, NM+1) = frac_bor(ij, ip, NM+1) &
                                               + frac_bor(ij, ip, im)*frac_phi(ij, ip, im)/sum(frac_phi(ij, ip, :))
                    frac_bor(ij, NP+1, im) = frac_bor(ij, NP+1, im) &
                                               + frac_bor(ij, ip, im)*frac_phi(ij, ip, im)/sum(frac_phi(ij, :, im))
                    frac_bor(ij, NP+1, NM+1) = frac_bor(ij, NP+1, NM+1) &
                                               + frac_bor(ij, ip, im)*frac_phi(ij, ip, im)/sum(frac_phi(ij, :, :))
                enddo
            enddo
        enddo

        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 1, 0), legend='Low Skilled - Good Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 1, 1), legend='Low Skilled - Bad Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 1, 2), legend='Low Skilled - Average')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 2, 0), legend='High Skilled - Good Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 2, 1), legend='High Skilled - Bad Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 2, 2), legend='High Skilled - Average')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 3, 0), legend='Good Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 3, 1), legend='Bad Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 3, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained Households')

        call plot(dble(ages), frac_phi(:, 1, 0)/(frac_phi(:, 1, 0) &
                                + frac_phi(:, 1, 1)), legend='Low Skilled - Good health')
        call plot(dble(ages), frac_phi(:, 1, 1)/(frac_phi(:, 1, 0) &
                                + frac_phi(:, 1, 1)), legend='Low Skilled - Bad health')
        call plot(dble(ages), frac_phi(:, 2, 0)/(frac_phi(:, 2, 0) &
                                + frac_phi(:, 2, 1)), legend='High Skilled - Good health')
        call plot(dble(ages), frac_phi(:, 2, 1)/(frac_phi(:, 2, 0) &
                                + frac_phi(:, 2, 1)), legend='High Skilled - Bad health')
        call plot(dble(ages), frac_phi(:, 1, 0) + frac_phi(:, 2, 0), legend='Average - Good health')
        call plot(dble(ages), frac_phi(:, 1, 1) + frac_phi(:, 2, 1), legend='Average - Bad health')
        call execplot(xlabel='Age j', ylabel='Frac. Good vs. Bad Health Households')

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, ip, im, is, ih

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do ip = 1, NP
                    do im = 0, NM
                        do is = 1, NS
                            do ih = 1, NH
                                if(phi(ij, ia, ip, im, is, ih) > 1d-8)iamax(ij) = ia
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
