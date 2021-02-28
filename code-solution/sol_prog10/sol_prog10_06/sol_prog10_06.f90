!##############################################################################
! PROGRAM Fertility
!
! ## The baseline life cycle model with endogenous fertility
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
include "sol_prog10_06m.f90"

program Fertility

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
        pen(JR:JJ) = 0.25d0*w*sum(eff)/dble(JR-1)

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi)
        eta = exp(eta)

        ! initialize probabilities for non-fertile period
        pf = 1d0
        call grid_Cons_grow(pf(1:JF, NF, 0), 0.02d0, 1.0d0, 0.2d0)
        pf(:, NF, NF) = 1d0 - pf(:, NF, 0)
        pf(:, 0, NF) = 1d0 - pf(:, 0, 0)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ia, in, is, is_max, ift, ift_max, ic
        real*8 :: x_in, wage
        real*8 :: c_temp(0:NC), aplus_temp(0:NC), u_temp(0:NC), uc_temp(0:NC), un_temp(0:NC)
        logical :: check = .false.

        ! get decision in the last period of life
        do ia = 0, NA
            chplus(JJ, ia, :, :, :) = 0
            aplus(JJ, ia, :, :, :) = 0d0
            c(JJ, ia, :, :, :) = (1d0+r)*a(ia) + pen(JJ)
            VC(JJ, ia, :, :, :) = valuefunc_c(0d0, c(JJ, ia, 0, 1, 0), 0, JJ, 0, 1, 0)
            VN(JJ, ia, :, :, :) = valuefunc_n(0d0, 0, JJ, 0, 1, 0)
        enddo

        ! interpolate individual RHS
        call interpolate(JJ)

        do ij = JJ-1, 1, -1

            ! check about how many is to iterate
            if(ij >= JR)then
                is_max = 1
                ift_max = 0
            else
                is_max = NS
                ift_max = NF
            endif

            do ia = 0, NA
                do in = 0, NN

                    ! determine decision for zero assets at retirement without pension
                    if(ij >= JR .and. ia == 0 .and. pen(ij) <= 1d-10)then
                        aplus(ij, ia, in, :, :) = 0d0
                        chplus(ij, ia, in, :, :) = 0
                        c(ij, ia, in, :, :) = 0d0
                        VC(ij, ia, in, :, :) = valuefunc_c(0d0, 0d0, 0, ij, 0, 1, 0)
                        VN(ij, ia, in, :, :) = valuefunc_n(0d0, 0, ij, 0, 1, 0)
                        cycle
                    endif

                    do is = 1, is_max
                        do ift = 0, ift_max

                            ! get initial guess for the individual choices
                            x_in = aplus(ij+1, ia, in, is, ift)

                            ! set up communication variables
                            ij_com = ij
                            ia_com = ia
                            in_com = in
                            is_com = is
                            ift_com = ift

                            if(ift == 0)then

                                ! set up communication variable
                                chplus_com = 0

                                ! solve the household problem using rootfinding
                                call fzero(x_in, foc, check)

                                ! check for borrowing constraint
                                if(x_in < 0d0)then
                                    x_in = 0d0
                                    wage = w*eff(ij)*eta(is)*(1d0-dble(in)**nu*lc)
                                    cons_com = (1d0+r)*a(ia) + wage + pen(ij)
                                endif

                                ! write screen output in case of a problem
                                if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, in, is, ift

                                aplus(ij, ia, in, is, ift) = x_in
                                c(ij, ia, in, is, ift) = cons_com
                                chplus(ij, ia, in, is, ift) = chplus_com
                                VC(ij, ia, in, is, ift) = valuefunc_c(x_in, cons_com, 0, ij, in, is, 0)
                                VN(ij, ia, in, is, ift) = valuefunc_n(x_in, 0, ij, in, is, 0)
                                V(ij, ia, in, is, ift) = VC(ij, ia, in, is, 0) + VN(ij, ia, in, is, 0)

                            else

                                do ic = 0, NC

                                    ! set up communication variable
                                    chplus_com = ic

                                    ! solve the household problem using rootfinding
                                    call fzero(x_in, foc, check)

                                    ! write screen output in case of a problem
                                    if(check)write(*,'(a, 5i4)')'ERROR IN ROOTFINDING : ', ij, ia, in, is, ift

                                    ! check for borrowing constraint
                                    if(x_in < 0d0)then
                                        x_in = 0d0
                                        wage = w*eff(ij)*eta(is)*(1d0-dble(in)**nu*lc)
                                        cons_com = (1d0+r)*a(ia) + wage + pen(ij)
                                    endif

                                    aplus_temp(ic) = x_in
                                    c_temp(ic) = cons_com
                                    uc_temp(ic) = valuefunc_c(x_in, cons_com, ic, ij, in, is, 1)
                                    un_temp(ic) = valuefunc_n(x_in, ic, ij, in, is, 1)
                                    u_temp(ic) = uc_temp(ic) + un_temp(ic)
                                enddo

                                ! choose the number of children that give more utility
                                if(u_temp(1) >= u_temp(0))then
                                    aplus(ij, ia, in, is, ift) = aplus_temp(1)
                                    c(ij, ia, in, is, ift) = c_temp(1)
                                    chplus(ij, ia, in, is, ift) = 1
                                    VC(ij, ia, in, is, ift) = uc_temp(1)
                                    VN(ij, ia, in, is, ift) = un_temp(1)
                                    V(ij, ia, in, is, ift) = u_temp(1)
                                else
                                    aplus(ij, ia, in, is, ift) = aplus_temp(0)
                                    c(ij, ia, in, is, ift) = c_temp(0)
                                    chplus(ij, ia, in, is, ift) = 0
                                    VC(ij, ia, in, is, ift) = uc_temp(0)
                                    VN(ij, ia, in, is, ift) = un_temp(0)
                                    V(ij, ia, in, is, ift) = u_temp(0)
                                endif
                            endif
                        enddo
                    enddo

                    ! copy decision in retirement age
                    if(ij >= JR)then
                        aplus(ij, ia, in, :, :) = aplus(ij, ia, in, 1, 0)
                        c(ij, ia, in, :, :) = c(ij, ia, in, 1, 0)
                        chplus(ij, ia, in, :, :) = chplus(ij, ia, in, 1, 0)
                        VC(ij, ia, in, :, :) = VC(ij, ia, in, 1, 0)
                        VN(ij, ia, in, :, :) = VN(ij, ia, in, 1, 0)
                        V(ij, ia, in, :, :) = V(ij, ia, in, 1, 0)
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
        integer :: ia, in, in_p, is, is_p, ift, ift_p
        real*8 :: chelp

        ! get RHS of first order condition
        do ia = 0, NA
            do in = 0, NN
                do is = 1, NS
                    do ift = 0, NF

                        ! calculate the RHS of the first order condition
                        RHS(ij, ia, in, is, ift) = 0d0
                        EVC(ij, ia, in, is, ift) = 0d0
                        EVN(ij, ia, in, is, ift) = 0d0
                        do in_p = 0, in
                            do is_p = 1, NS
                                do ift_p = 0, NF
                                    chelp = max(c(ij, ia, in_p, is_p, ift_p), 1d-10)
                                    RHS(ij, ia, in, is, ift) = RHS(ij, ia, in, is, ift) + &
                                      pi(is, is_p)*pf(ij, ift, ift_p)*binomialPDF(in_p, in, p)* &
                                      margu(chelp)
                                    EVC(ij, ia, in, is, ift) = EVC(ij, ia, in, is, ift) + &
                                      pi(is, is_p)*pf(ij, ift, ift_p)*binomialPDF(in_p, in, p)* &
                                      VC(ij, ia, in_p, is_p, ift_p)
                                    EVN(ij, ia, in, is, ift) = EVN(ij, ia, in, is, ift) + &
                                      pi(is, is_p)*pf(ij, ift, ift_p)*binomialPDF(in_p, in, p)* &
                                      VN(ij, ia, in_p, is_p, ift_p)
                                enddo
                            enddo
                        enddo
                        RHS(ij, ia, in, is, ift) = ((1d0+r)*beta*psi(ij)*RHS(ij, ia, in, is, ift))**(-gamma)
                        EVC(ij, ia, in, is, ift) = (egam*EVC(ij, ia, in, is, ift))**(1d0/egam)
                        EVN(ij, ia, in, is, ift) = (echi*EVN(ij, ia, in, is, ift))**(1d0/echi)
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! determines the invariant distribution of households
    subroutine get_distribution()

        implicit none
        integer :: ij, ia, in, is, ift, ift_p, is_p, in_p, ic, in_h, ial, iar
        real*8 :: varphi

        ! set distribution to zero
        phi(:, :, :, :, :) = 0d0

        ! get initial distribution in age 1
        do ift_p = 0, NF
            phi(1, 0, 0, is_initial, ift_p) = pf(1, NF, ift_p)
        enddo

        ! successively compute distribution over ages
        do ij = 2, JJ

            ! iterate over yesterdays gridpoints
            do ia = 0, NA
                do in = 0, NN
                    do is = 1, NS
                        do ift = 0, NF

                            ! interpolate yesterday's savings decision
                            call linint_Grow(aplus(ij-1, ia, in, is, ift), a_l, a_u, a_grow, NA, ial, iar, varphi)

                            ! restrict values to grid just in case
                            ial = min(ial, NA)
                            iar = min(iar, NA)
                            varphi = min(varphi, 1d0)

                            ! redistribute households
                            ic = int(chplus(ij-1, ia, in, is, ift))
                            in_h = in + ic
                            if(in_h > NN)in_h = NN

                            ! redistribute households
                            do is_p = 1, NS
                                do ift_p = 0, NF
                                    do in_p = 0, in_h
                                        phi(ij, ial, in_p, is_p, ift_p) = phi(ij, ial, in_p, is_p, ift_p) &
                                            + pf(ij, ift, ift_p)*pi(is, is_p)*binomialPDF(in_p, in_h, p) &
                                            *varphi*phi(ij-1, ia, in, is, ift)
                                        phi(ij, iar, in_p, is_p, ift_p) = phi(ij, iar, in_p, is_p, ift_p) &
                                            + pf(ij, ift, ift_p)*pi(is, is_p)*binomialPDF(in_p, in_h, p) &
                                            *(1d0-varphi)*phi(ij-1, ia, in, is, ift)
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
        integer :: ij, ia, in, is, ift

        ! calculate cohort averages
        c_coh(:) = 0d0
        y_coh(:) = 0d0
        a_coh(:) = 0d0
        v_coh(:) = 0d0
        ch_coh(:) = 0d0
        n_ch_coh(:) = 0d0
        age_ch_coh(:) = 0d0
        sum_phi(:) = 0d0

        do ij = 1, JJ
            do ia = 0, NA
                do in = 0, NN
                    do is = 1, NS
                        do ift = 0, NF
                            c_coh(ij) = c_coh(ij) + c(ij, ia, in, is, ift)*phi(ij, ia, in, is, ift)
                            y_coh(ij) = y_coh(ij) + w*eff(ij)*eta(is)**(1d0-in**nu*lc)*phi(ij, ia, in, is, ift)
                            a_coh(ij) = a_coh(ij) + a(ia)*phi(ij, ia, in, is, ift)
                            if(chplus(ij, ia, in, is, ift) == 1)then
                                n_ch_coh(ij) = n_ch_coh(ij) + phi(ij, ia, in, is, ift)
                                age_ch_coh(in) = age_ch_coh(in) + dble(ij+1)*phi(ij, ia, in, is, ift)
                                sum_phi(in) = sum_phi(in) + phi(ij, ia, in, is, ift)
                            endif
                            ch_coh(ij) = ch_coh(ij) + dble(in)*phi(ij, ia, in, is, ift)
                            v_coh(ij) = v_coh(ij) + V(ij, ia, in, is, ift)*phi(ij, ia, in, is, ift)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        ! compute average age at which the n-th child is born
        do in = 1, NN
            if(age_ch_coh(in) > 0)age_ch_coh(in) = age_ch_coh(in)/sum_phi(in)
        enddo

        ! calculate cohort specific coeffcients of variation
        cv_c(:) = 0d0
        cv_y(:) = 0d0
        do ij = 1, JJ
            do ia = 0, NA
                do in = 0, NN
                    do is = 1, NS
                        do ift = 0, NF
                            cv_c(ij) = cv_c(ij) + c(ij, ia, in, is, ift)**2*phi(ij, ia, in, is, ift)
                            cv_y(ij) = cv_y(ij) + (eff(ij)*eta(is))**2*phi(ij, ia, in, is, ift)
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
        integer :: ift, ij, ia, in, is, iamax(JJ), inmax(JJ), ages(JJ)
        real*8 :: frac_bor(JJ)

        ! check for the maximium grid points used
        call check_grid_a(iamax)
        call check_grid_n(inmax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        write(21,'(a,a)')' IJ      CONS  EARNINGS    INCOME      PENS    ASSETS  N_BIRTHS    C_HOME', &
            '     CV(C)     CV(L)     VALUE     IAMAX     INMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,2i10)')ij, c_coh(ij), y_coh(ij), &
                y_coh(ij)+r*a_coh(ij), pen(ij), a_coh(ij), sum(n_ch_coh(1:ij)), ch_coh(ij), &
                cv_c(ij), cv_y(ij), v_coh(ij), iamax(ij), inmax(ij)
        enddo
        write(21, '(a/)') '--------------------------------------------------------------------'
        write(*, '(a, f10.2)')'Number of children at Age 30:', sum(n_ch_coh(1:10))
        write(*, '(a, f10.2)')'Number of children at Age 45:', sum(n_ch_coh(1:25))
        write(*, '(a, f10.2)')'Age at birth of first child: ', 20 + age_ch_coh(0)
        write(*, '(a, f10.2)')'Age at birth of second child:', 20 + age_ch_coh(1)
        write(*, '(a, f10.2/)')'Children at home at Age 37:  ', ch_coh(37-20)

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption')
        call plot(dble(ages), y_coh+pen, legend='Earnings')
        call execplot(xlabel='Age j', ylabel='Mean')

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

        call plot(dble(ages), cv_c, legend='Consumption')
        call plot(dble(ages), cv_y, legend='Earnings')
        call execplot(xlabel='Age j', ylabel='Coefficient of Variation')

        call plot(dble(ages), ch_coh, legend='Children at Home')
        call plot(dble(ages), n_ch_coh, legend='Newborn Children')
        call execplot(xlabel='Age j', ylabel='Mean')

        ! calculate fraction of borrowing constrained households
        frac_bor = 0d0
        do ij = 1, JJ-1
            do ift = 0, NF
                do ia = 0, NA
                    do in = 0, NN
                        do is = 1, NS
                            if(aplus(ij, ia, in, is, ift) < 1d-6)then
                                frac_bor(ij) = frac_bor(ij) + phi(ij, ia, in, is, ift)
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
        frac_bor(JJ) = 1d0

        call plot(dble(ages(1:JJ-1)), frac_bor(1:JJ-1))
        call execplot(xlabel='Age j', ylabel='Frac. Borrowing Constrained Households')

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, in, is, ift

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do in = 0, NN
                    do is = 1, NS
                        do ift = 0, 1
                            if(phi(ij, ia, in, is, ift) > 1d-8)iamax(ij) = ia
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum child gridpoint used
    subroutine check_grid_n(inmax)

        implicit none
        integer :: inmax(JJ), ij, ia, in, is, ift

        inmax = 0
        do ij = 1, JJ

            ! check for the maximum child grid point used at a certain age
            do in = 0, NN
                do ia = 0, NA
                    do is = 1, NS
                        do ift = 0, 1
                            if(phi(ij, ia, in, is, ift) > 1d-8)inmax(ij) = in
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine

end program
