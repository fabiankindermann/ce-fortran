!##############################################################################
! PROGRAM LaborHealth
!
! ## The baseline life cycle model with labor supply depending on health risk
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog10_08m.f90"

program LaborHealth

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

        ! health stauts affects survival probabilities
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
        pen(JR:JJ) = 0.5d0*sum(eff)/dble(JR-1)*0.33d0

        ! initialize health shock
        dist_m(0) = 0.95d0
        dist_m(1) = 0.05d0

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! probability to have bad health when current health is good
        call grid_Cons_Equi(pi_m(:, 0, 1), 0.1d0, 0.4d0)
        pi_m(:, 0, 0) = 1d0 - pi_m(:, 0, 1)

        ! probability to have bad health when current health is bad
        call grid_Cons_Equi(pi_m(:, 1, 1), 0.6d0, 0.9d0)
        pi_m(:, 1, 0) = 1d0 - pi_m(:, 1, 1)

        ! initialize impact of shock on productivity
        varrho_m(0) = 1d0
        varrho_m(1) = exp(-varrho)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, im, is, is_max
        real*8 :: x_in, wage, available
        logical :: check

        ! get decision in the last period of life
        do ia = 0, NA
            do im = 0, NM
                aplus(JJ, ia, im, :) = 0d0
                c(JJ, ia, im, :) = (1d0+r)*a(ia) + pen(JJ)
                l(JJ, ia, im, :) = 0d0
                V(JJ, ia, im, :) = valuefunc(0d0, c(JJ, ia, im, 1), l(JJ, ia, im, 1), JJ, im, 1)
            enddo
        enddo

        ! interpolate individual RHS
        call interpolate(JJ)

        do ij = JJ-1, 1, -1

            ! check about how many is to iterate
            if(ij >= JR)then
                is_max = 1
            else
                is_max = NS
            endif

            do ia = 0, NA
                do im = 0, NM

                    ! determine decision for zero assets at retirement without pension
                    if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                        aplus(ij, ia, im, :) = 0d0
                        c(ij, ia, im, :) = 0d0
                        l(ij, ia, im, :) = 0d0
                        V(ij, ia, im, :) = valuefunc(0d0, 0d0, 0d0, ij, im, 1)
                        cycle
                    endif

                    do is = 1, is_max

                        ! get initial guess for the individual choices
                        x_in = aplus(ij+1, ia, im, is)

                        ! set up communication variables
                        ij_com = ij
                        ia_com = ia
                        im_com = im
                        is_com = is

                        ! solve the household problem using rootfinding
                        call fzero(x_in, foc, check)

                        ! write screen output in case of a problem
                        if(check)write(*,'(a, 4i4)')'ERROR IN ROOTFINDING : ', ij, ia, im, is

                        ! check for borrowing constraint
                        if(x_in < 0d0)then
                            x_in = 0d0
                            wage = w*eff(ij)*eta(is)*varrho_m(im)
                            available = (1d0+r)*a(ia) + pen(ij)
                            if(ij < JR)then
                                lab_com = min(max(nu*(1d0-phi_l*dble(im)) - (1d0-nu)*available/wage, 1d-10) &
                                              , 1d0-phi_l*dble(im)-1d-10)
                            else
                                lab_com = 0d0
                            endif
                            cons_com = max((available + wage*lab_com), 1d-10)
                        endif

                        ! copy decisions
                        aplus(ij, ia, im, is) = x_in
                        c(ij, ia, im, is) = cons_com
                        l(ij, ia, im, is) = lab_com
                        V(ij, ia, im, is) = valuefunc(x_in, cons_com, lab_com, ij, im, is)
                    enddo

                    ! copy decision in retirement age
                    if(ij >= JR)then
                        aplus(ij, ia, im, :) = aplus(ij, ia, im, 1)
                        c(ij, ia, im, :) = c(ij, ia, im, 1)
                        l(ij, ia, im, :) = l(ij, ia, im, 1)
                        V(ij, ia, im, :) = V(ij, ia, im, 1)
                    endif
                enddo
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
        integer :: ia, im, im_p, is, is_p
        real*8 :: chelp, lhelp

        do ia = 0, NA
            do im = 0, NM
                do is = 1, NS

                    ! calculate the RHS of the first order condition
                    RHS(ij, ia, im, is) = 0d0
                    EV(ij, ia, im, is) = 0d0
                    do im_p = 0, NM
                        do is_p = 1, NS
                            chelp = max(c(ij, ia, im_p, is_p), 1d-10)
                            lhelp = max(l(ij, ia, im_p, is_p), 1d-10)
                            RHS(ij, ia, im, is) = RHS(ij, ia, im, is) + &
                                                  pi_m(ij, im, im_p)*pi(is, is_p)*margu(chelp, lhelp, im_p)
                            EV(ij, ia, im, is) = EV(ij, ia, im, is) + &
                                                 pi_m(ij, im, im_p)*pi(is, is_p)*V(ij, ia, im_p, is_p)
                        enddo
                    enddo
                    RHS(ij, ia, im, is) = ((1d0+r)*beta*psi(ij, im)*RHS(ij, ia, im, is))**(-gamma)
                    EV(ij, ia, im, is) = (egam*EV(ij, ia, im, is))**(1d0/egam)
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, im, im_p, is, is_p, ial, iar
        real*8 :: varphi

        ! set distribution to zero
        phi(:, :, :, :) = 0d0

        ! get initial distribution in age 1
        do im = 0, NM
            phi(1, 0, im, is_initial) = dist_m(im)
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do im = 0, NM
                    do is = 1, NS

                        ! interpolate yesterday's savings decision
                        call linint_Grow(aplus(ij-1, ia, im, is), a_l, a_u, a_grow, NA, ial, iar, varphi)

                        ! restrict values to grid just in case
                        ial = min(ial, NA)
                        iar = min(iar, NA)
                        varphi = min(varphi, 1d0)

                        ! redistribute households
                        do im_p = 0, NM
                            do is_p = 1, NS
                                phi(ij, ial, im_p, is_p) = phi(ij, ial, im_p, is_p) + &
                                    pi_m(ij, im, im_p)*pi(is, is_p)*varphi*phi(ij-1, ia, im, is)
                                phi(ij, iar, im_p, is_p) = phi(ij, iar, im_p, is_p) + &
                                    pi_m(ij, im, im_p)*pi(is, is_p)*(1d0-varphi)*phi(ij-1, ia, im, is)
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
        integer :: ij, ia, im, is

        ! calculate fraction of good vs. bad health households
        do ij = 1, JJ
            do im = 0, NM
                frac_phi(ij, im) = sum(phi(ij, :, im, :))
            enddo
        enddo

        ! calculate cohort averages
        c_coh(:, :) = 0d0
        y_coh(:, :) = 0d0
        a_coh(:, :) = 0d0
        l_coh(:, :) = 0d0
        v_coh(:, :) = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do im = 0, NM
                    do is = 1, NS
                        c_coh(ij, im) = c_coh(ij, im) + c(ij, ia, im, is) &
                                        *phi(ij, ia, im, is)/frac_phi(ij, im)
                        y_coh(ij, im) = y_coh(ij, im) + eff(ij)*eta(is)*varrho_m(im) &
                                        *phi(ij, ia, im, is)/frac_phi(ij, im)
                        a_coh(ij, im) = a_coh(ij, im) + a(ia) &
                                        *phi(ij, ia, im, is)/frac_phi(ij, im)
                        l_coh(ij, im) = l_coh(ij, im) + l(ij, ia, im, is) &
                                        *phi(ij, ia, im, is)/frac_phi(ij, im)
                        v_coh(ij, im) = v_coh(ij, im) + V(ij, ia, im, is) &
                                        *phi(ij, ia, im, is)/frac_phi(ij, im)
                    enddo
                enddo
            enddo
        enddo

        ! recover unconditional cohort averages
        do ij = 1, JJ
            do im = 0, NM
                c_coh(ij, NM+1) = c_coh(ij, NM+1) + c_coh(ij, im)*frac_phi(ij, im)
                y_coh(ij, NM+1) = y_coh(ij, NM+1) + y_coh(ij, im)*frac_phi(ij, im)
                a_coh(ij, NM+1) = a_coh(ij, NM+1) + a_coh(ij, im)*frac_phi(ij, im)
                l_coh(ij, NM+1) = l_coh(ij, NM+1) + l_coh(ij, im)*frac_phi(ij, im)
                v_coh(ij, NM+1) = v_coh(ij, NM+1) + v_coh(ij, im)*frac_phi(ij, im)
            enddo
        enddo

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ia, im, is, iamax(JJ), ages(JJ)
        real*8 :: frac_bor(JJ, 0:NM+1)

        ! check for the maximium grid point used
        call check_grid(iamax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        ! set up output file
        do im = 0, NM+1

            if(im == 0)then
                write(21,'(a/)') 'GOOD HEALTH'
            elseif(im == 1)then
                write(21,'(a/)') 'BAD HEALTH'
            else
                write(21,'(a/)') 'UNCONDITIONAL COHORT AVERAGE'
            endif

            write(21,'(a,a)') ' IJ      CONS  EARNINGS    INCOME      PENS    ASSETS', &
                '     VALUE     IAMAX'
            do ij = 1, JJ
                write(21,'(i3,6f10.3,i10)') ij, c_coh(ij, im), y_coh(ij, im), &
                    y_coh(ij, im)+r*a_coh(ij, im), pen(ij), a_coh(ij, im), &
                    v_coh(ij, im), iamax(ij)
            enddo
            write(21,'(a/)') '--------------------------------------------------------------------'

        enddo

        ! plot output
        call plot(dble(ages), c_coh(:, 0), legend='Good Health')
        call plot(dble(ages), c_coh(:, 1), legend='Bad Health')
        call plot(dble(ages), c_coh(:, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Consumption')

        call plot(dble(ages), l_coh(:, 0), legend='Good health')
        call plot(dble(ages), l_coh(:, 1), legend='Bad health')
        call plot(dble(ages), l_coh(:, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Working time')

        call plot(dble(ages), y_coh(:, 0)+pen, legend='Good Health')
        call plot(dble(ages), y_coh(:, 1)+pen, legend='Bad Health')
        call plot(dble(ages), y_coh(:, 2)+pen, legend='Average')
        call execplot(xlabel='Age j', ylabel='Earnings')

        call plot(dble(ages), a_coh(:, 0), legend='Good Health')
        call plot(dble(ages), a_coh(:, 1), legend='Bad Health')
        call plot(dble(ages), a_coh(:, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Assets')

        ! calculate fraction of borrowing constrained households
        frac_bor = 0d0
        do ij = 1, JJ-1
            do ia = 0, NA
                do im = 0, NM
                    do is = 1, NS
                        if(aplus(ij, ia, im, is) < 1d-6)then
                            frac_bor(ij, im) = frac_bor(ij, im) + phi(ij, ia, im, is)/frac_phi(ij, im)
                            frac_bor(ij, NM+1) = frac_bor(ij, NM+1) + phi(ij, ia, im, is)
                        endif
                    enddo
                enddo
            enddo
        enddo
        frac_bor(JJ, :) = 1d0

        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 0), legend='Good Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 1), legend='Bad Health')
        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1, 2), legend='Average')
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained Households')

        call plot(dble(ages), frac_phi(:, 0), legend='Good Health')
        call plot(dble(ages), frac_phi(:, 1), legend='Bad Health')
        call execplot(xlabel='Age j', ylabel='Frac. Good vs. Bad Health Households')

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, im, is

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do im = 0, NM
                    do is = 1, NS
                        if(phi(ij, ia, im, is) > 1d-8)iamax(ij) = ia
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
