!##############################################################################
! PROGRAM PortfolioChoiceAnnuities
!
! ## Portfolio choice with the option to annuitize wealth upon retirement
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
include "prog10_05m.f90"

program PortfolioChoiceAnnuities

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
        real*8 :: temp(NSR, 2), psi_cum
        integer :: ij

        ! wage rate for effective labor
        w = 1d0

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
        pen(JR:JJ) = kappa*w*eff(JR-1)

        ! discretize zeta shocks
        call normal_discrete(zeta, dist_zeta, 0d0, sigma_zeta)
        zeta = exp(zeta)

        ! discretize eps-vtheta shocks
        call normal_discrete((/NS, NR/), temp, dist_epsvtheta, (/0d0, 0d0/), &
                                                  (/sigma_eps, sigma_vtheta/), rho)
        eps(:) = exp(temp(:, 1))
        vtheta(:)  = temp(:, 2)

        ! initialize asset grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize annuity value grid
        call grid_Cons_Grow(ar, ar_l, ar_u, ar_grow)

        ! calculate the price of an annuity
        p_a = 0d0
        psi_cum = 1d0
        do ij = JR, JJ
            psi_cum = psi_cum*psi(ij)
            p_a = p_a + psi_cum/(1d0+r_f)**(ij-(JR-1))
        enddo
        p_a = (1d0+xi)*p_a

        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = min(w*minval(eff(1:JR-1))*minval(eps(:))*zeta(1), pen(JR))
        X_u = (1d0 + r_f + mu_r + maxval(vtheta(:)))*a_u + &
                    max(w*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW), ar_u/p_a)
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ix, ia, ir

        ! get decisions in last period of life
        omega_plus(JJ, :, :) = 0d0
        V_ar(JJ, :, :) = 0d0
        do ix = 0, NX
            a_plus(JJ, ix, :) = 0d0
            c(JJ, ix, :) = X(ix)
            V(JJ, ix, :) = valuefunc(0d0, c(JJ, ix, 0), JJ, 0)
        enddo

        ! solve household problem up to date of retirement
        do ij = JJ-1, JR, -1

            do ia = 1, NA
                do ir = 0, NAR
                    call solve_portfolio(ij, ia, ir)
                enddo
            enddo
            omega_plus(ij, 0, :) = omega_plus(ij, 1, :)

            call interpolate(ij)

            do ix = 0, NX
                do ir = 0, NAR
                    call solve_consumption(ij, ix, ir)
                enddo
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

        ! solve portfolio choice problem directly before retirement
        do ia = 1, NA
            do ir = 0, NAR
                call solve_portfolio(ij, ia, ir)
            enddo
        enddo
        omega_plus(JR-1, 0, :) = omega_plus(JR-1, 1, :)

        ! get optimal fraction to annuitize
        do ia = 0, NA
            call solve_annuitization(ia)
        enddo

        ! now interpolate using information about annuity choice
        call interpolate(JR-1)

        ! solve the consumption savings problem
        do ix = 0, NX
            call solve_consumption(JR-1, ix, 0)
        enddo

        write(*,'(a,i3,a)')'Age: ',JR-1,' DONE!'

        ! now solve remainder of the life cycle
        do ij = JR-2, 1, -1

            do ia = 1, NA
                call solve_portfolio(ij, ia, 0)
            enddo
            omega_plus(ij, 0, 0) = omega_plus(ij, 1, 0)

            call interpolate(ij)

            do ix = 0, NX
                call solve_consumption(ij, ix, 0)
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia, ir)

        implicit none
        integer, intent(in) :: ij, ia, ir
        real*8 :: x_in, port0, port1, tolerance
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        ir_com = ir

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

        ! use intermediate value theorem
        if(port0*port1 > 0d0)then
            if(abs(port0) > abs(port1))then
                omega_plus(ij, ia, ir) = 1d0
            else
                omega_plus(ij, ia, ir) = 0d0
            endif
            return
        else

            ! get order of magnitude of foc
            tolerance = 1d-5*abs(port0-port1)
            tolerance = min(tolerance, 1d-8)
            call settol_root(tolerance)

            ! get best guess for the root of foc_port
            x_in = -port0/(port1-port0)
            check = .false.

            ! solve the household problem using rootfinding
            call fzero(x_in, foc_port, check)

            ! write screen output in case of a problem
            if(check)write(*,'(a, 3i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia, ir

            omega_plus(ij, ia, ir) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ix, ir)

        implicit none
        integer, intent(in) :: ij, ix, ir
        real*8 :: x_in
        logical :: check

        ! determine decision for zero cash-on-hand
        if(X(ix) < 1d-10)then
            a_plus(ij, ix, ir) = 0d0
            c(ij, ix, ir) = 0d0
            V(ij, ix, ir) = valuefunc(0d0, 0d0, ij, ir)
            if(ij >= JR)V_ar(ij, ix, ir) = dV_dar(0d0, ij, ir)
            return
        endif

        ! set up communication variables
        ij_com = ij
        ix_com = ix
        ir_com = ir

        ! get best initial guess from future period
        x_in = a_plus(ij+1, ix, ir)
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 3i4)')'ERROR IN ROOTFINDING CONS : ', ij, ix, ir

        ! check for borrowing constraint
        if(x_in < 0d0)then
            x_in = 0d0
            cons_com = X(ix)
        endif

        ! copy decisions
        a_plus(ij, ix, ir) = x_in
        c(ij, ix, ir) = cons_com
        V(ij, ix, ir) = valuefunc(x_in, cons_com, ij, ir)
        if(ij >= JR)V_ar(ij, ix, ir) = dV_dar(x_in, ij, ir)

    end subroutine


    ! solve the household's annuitization decision
    subroutine solve_annuitization(ia)

        implicit none
        integer, intent(in) :: ia
        real*8 :: x_in, annu0, annu1
        logical :: check

        ! set up communication variables
        ia_com = ia

        ! check for corner solutions
        annu0 = foc_annu(0d0)
        annu1 = foc_annu(a(ia))

        ! use intermediate value theorem
        if(annu0*annu1 > 0d0)then
            if(abs(annu0) > abs(annu1))then
                omegar_plus(ia) = 1d0
            else
                omegar_plus(ia) = 0d0
            endif
            return
        else

            ! get best guess for the root of foc_port
            if(ia > 0)then
                x_in = omegar_plus(ia-1)
            else
                x_in = 0d0
            endif
            check = .false.

            ! solve the household problem using rootfinding
            call fzero(x_in, foc_annu, check)

            ! write screen output in case of a problem
            if(check)write(*,'(a, i4)')'ERROR IN ROOTFINDING ANNU : ', ia

            omegar_plus(ia) = x_in
        endif

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, iw, isr, ir
        real*8 :: X_p, c_p, varphi_X, dist, EV, R_port
        real*8 :: ar_p, al_p, varphi_a, varphi_r, omega_p, Var_p
        integer :: ixl, ixr, ial, iar, irl, irr

        RHS(ij, :, :) = 0d0
        Q(ij, :, :) = 0d0

        ! case agent is retired
        if(ij >= JR)then
            do ia = 0, NA
                do ir = 0, NAR
                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, ir)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia) + pen(ij+1) + ar(ir)/p_a
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! get future consumption value
                        c_p = varphi_X      *c(ij+1, ixl, ir) + &
                              (1d0-varphi_X)*c(ij+1, ixr, ir)
                        c_p = max(c_p, 1d-10)

                        ! get tomorrow's value function
                        EV = varphi_X      *(egam*V(ij+1, ixl, ir))**(1d0/egam) + &
                             (1d0-varphi_X)*(egam*V(ij+1, ixr, ir))**(1d0/egam)

                        ! get RHS of foc and Q
                        RHS(ij, ia, ir) = RHS(ij, ia, ir) + dist*R_port*margu(c_p)
                        Q(ij, ia, ir)   = Q(ij, ia, ir)   + dist*EV**egam/egam
                    enddo

                    RHS(ij, ia, ir) = (beta*psi(ij+1)*RHS(ij, ia, ir))**(-gamma)
                    Q(ij, ia, ir)   = (egam*Q(ij, ia, ir))**(1d0/egam)
                enddo
            enddo

        ! special case of annuity choice
        elseif(ij == JR-1)then

            S(:) = 0d0

            do ia = 0, NA

                ! determine future liquid wealth and amount of annuity
                al_p = (1d0-omegar_plus(ia))*a(ia)
                ar_p = omegar_plus(ia)*a(ia)

                ! derive interpolation weights
                call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                call linint_Grow(ar_p, ar_l, ar_u, ar_grow, NAR, irl, irr, varphi_r)

                ! get optimal investment strategy
                if(varphi_a <= varphi_r)then
                    omega_p = varphi_a           *omega_plus(JR-1, ial, irl) + &
                              (varphi_r-varphi_a)*omega_plus(JR-1, iar, irl) + &
                              (1d0-varphi_r)     *omega_plus(JR-1, iar, irr)
                else
                    omega_p = varphi_r           *omega_plus(JR-1, ial, irl) + &
                              (varphi_a-varphi_r)*omega_plus(JR-1, ial, irr) + &
                              (1d0-varphi_a)     *omega_plus(JR-1, iar, irr)
                endif

                do isr = 1, NSR

                    ! get future cash-on-hand and interpolate
                    R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                    X_p = R_port*al_p+ pen(ij+1) + ar_p/p_a
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                    ! get distributional weight
                    dist = dist_epsvtheta(isr)

                    ! get future interpolated values
                    if(varphi_X <= varphi_r)then
                        c_p   = varphi_X           *c(ij+1, ixl, irl) + &
                                (varphi_r-varphi_X)*c(ij+1, ixr, irl) + &
                                (1d0-varphi_r)     *c(ij+1, ixr, irr)
                        Var_p = varphi_X           *V_ar(ij+1, ixl, irl) + &
                                (varphi_r-varphi_X)*V_ar(ij+1, ixr, irl) + &
                                (1d0-varphi_r)     *V_ar(ij+1, ixr, irr)
                        EV    = varphi_X           *(egam*V(ij+1, ixl, irl))**(1d0/egam) + &
                                (varphi_r-varphi_X)*(egam*V(ij+1, ixr, irl))**(1d0/egam) + &
                                (1d0-varphi_r)     *(egam*V(ij+1, ixr, irr))**(1d0/egam)
                    else
                        c_p   = varphi_r           *c(ij+1, ixl, irl) + &
                                (varphi_X-varphi_r)*c(ij+1, ixl, irr) + &
                                (1d0-varphi_X)     *c(ij+1, ixr, irr)
                        Var_p = varphi_r           *V_ar(ij+1, ixl, irl) + &
                                (varphi_X-varphi_r)*V_ar(ij+1, ixl, irr) + &
                                (1d0-varphi_X)     *V_ar(ij+1, ixr, irr)
                        EV    = varphi_r           *(egam*V(ij+1, ixl, irl))**(1d0/egam) + &
                                (varphi_X-varphi_r)*(egam*V(ij+1, ixl, irr))**(1d0/egam) + &
                                (1d0-varphi_X)     *(egam*V(ij+1, ixr, irr))**(1d0/egam)
                    endif
                    c_p = max(c_p, 1d-10)
                    Var_p = max(Var_p, 1d-10)**(-1d0/gamma)

                    ! get RHS of foc and S
                    RHS(ij, ia, 0) = RHS(ij, ia, 0) + dist*(margu(c_p)/p_a + Var_p)
                    S(ia)          = S(ia) + dist*EV**egam/egam
                enddo

                RHS(ij, ia, 0) = (beta*psi(ij+1)*RHS(ij, ia, 0))**(-gamma)
                S(ia)          = (egam*S(ia))**(1d0/egam)
            enddo

        ! agent is working
        else
            do ia = 0, NA
                do iw = 1, NW
                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, 0)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia)/eps(isr) + w*eff(ij+1)*zeta(iw)
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get distributional weight
                        dist = dist_zeta(iw)*dist_epsvtheta(isr)

                        ! get future consumption value
                        c_p = varphi_X      *c(ij+1, ixl, 0) + &
                              (1d0-varphi_X)*c(ij+1, ixr, 0)
                        c_p = max(c_p, 1d-10)

                        ! get tomorrow's value function
                        EV = varphi_X      *(egam*V(ij+1, ixl, 0))**(1d0/egam) + &
                             (1d0-varphi_X)*(egam*V(ij+1, ixr, 0))**(1d0/egam)

                        ! get RHS of foc and Q
                        RHS(ij, ia, 0) = RHS(ij, ia, 0) + dist*R_port*margu(eps(isr)*c_p)
                        Q(ij, ia, 0)   = Q(ij, ia, 0)   + dist*(eps(isr)*EV)**egam/egam
                    enddo
                enddo

                RHS(ij, ia, 0) = (beta*psi(ij+1)*RHS(ij, ia, 0))**(-gamma)
                Q(ij, ia, 0)   = (egam*Q(ij, ia, 0))**(1d0/egam)
            enddo
        endif

    end subroutine


    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij

        ! set distributions to zero
        phi_X(:, :, :) = 0d0
        phi_a(:, :, :) = 0d0

        do ij = 1, JJ

            ! get distribution on cash-on-hand grid
            call get_distribution_X(ij)

            ! get distribution on asset grid
            call get_distribution_a(ij)
        enddo

    end subroutine


    ! to calculate distribution on cash-on-hand grid
    subroutine get_distribution_X(ij)

        implicit none
        integer :: ij, ia, iw, ir, isr, ixl, ixr
        real*8 :: varphi, X_p, R_port, dist

        if(ij == 1)then

            ! get initial distribution at age 1 of cash-on-hand
            do iw = 1, NW

                ! get initial cash-on-hand
                X_p = w*eff(1)*zeta(iw)

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                ! get distributional weight
                dist = dist_zeta(iw)

                ! initialize the distribution
                phi_X(1, ixl, 0) = phi_X(1, ixl, 0) + dist*varphi
                phi_X(1, ixr, 0) = phi_X(1, ixr, 0) + dist*(1d0-varphi)
            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays asset distribution
            do ia = 0, NA

                ! iterate over current shocks
                do iw = 1, NW
                    do isr = 1, NSR

                        ! get today's cash-on-hand
                        R_port = 1d0 + r_f + omega_plus(ij-1, ia, 0)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia)/eps(isr) + w*eff(ij)*zeta(iw)

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                        ! get distributional weight
                        dist = dist_zeta(iw)*dist_epsvtheta(isr)

                        phi_X(ij, ixl, 0) = phi_X(ij, ixl, 0) + dist*varphi*phi_a(ij-1, ia, 0)
                        phi_X(ij, ixr, 0) = phi_X(ij, ixr, 0) + dist*(1d0-varphi)*phi_a(ij-1, ia, 0)
                    enddo
                enddo
            enddo

        else

            do ia = 0, NA
                do ir = 0, NAR

                    ! iterate over current shocks
                    do isr = 1, NSR

                        ! get today's cash-on-hand
                        R_port = 1d0 + r_f + omega_plus(ij-1, ia, ir)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia) + pen(ij) + ar(ir)/p_a

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        phi_X(ij, ixl, ir) = phi_X(ij, ixl, ir) + dist*varphi*phi_a(ij-1, ia, ir)
                        phi_X(ij, ixr, ir) = phi_X(ij, ixr, ir) + dist*(1d0-varphi)*phi_a(ij-1, ia, ir)
                    enddo
                enddo
            enddo
        endif

    end subroutine


    ! to calculate the end of period asset distribution
    subroutine get_distribution_a(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ix, ir, ial, iar, irl, irr
        real*8 :: varphi_a, varphi_r, ar_p, al_p, omega_annu

        if(ij < JR-1)then

            ! iterate over todays cash on hand
            do ix = 0, NX

                ! interpolate asset decision
                call linint_Grow(a_plus(ij, ix, 0), a_l, a_u, a_grow, NA, ial, iar, varphi_a)

                ! restrict values to grid just in case
                ial = min(ial, NA)
                iar = min(iar, NA)
                varphi_a = min(varphi_a, 1d0)

                ! get end of period asset distribution
                phi_a(ij, ial, 0) = phi_a(ij, ial, 0) + varphi_a*phi_X(ij, ix, 0)
                phi_a(ij, iar, 0) = phi_a(ij, iar, 0) + (1d0-varphi_a)*phi_X(ij, ix, 0)
            enddo

        elseif(ij == JR-1)then

            ! iterate over todays cash on hand
            do ix = 0, NX

                ! determine future annuitization share through interpolation
                call linint_Grow(a_plus(ij, ix, 0), a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                omega_annu = varphi_a*omegar_plus(ial) + (1d0-varphi_a)*omegar_plus(iar)

                ! determine future liquid wealth and amount of annuity
                al_p = (1d0-omega_annu)*a_plus(ij, ix, 0)
                ar_p = omega_annu*a_plus(ij, ix, 0)

                ! derive interpolation weights
                call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                call linint_Grow(ar_p, ar_l, ar_u, ar_grow, NAR, irl, irr, varphi_r)

                ! restrict values to grid just in case
                ial = min(ial, NA)
                iar = min(iar, NA)
                varphi_a = min(varphi_a, 1d0)

                ! get end of period asset distribution
                phi_a(ij, ial, irl) = phi_a(ij, ial, irl) + varphi_a*varphi_r*phi_X(ij, ix, 0)
                phi_a(ij, iar, irl) = phi_a(ij, iar, irl) + (1d0-varphi_a)*varphi_r*phi_X(ij, ix, 0)
                phi_a(ij, ial, irr) = phi_a(ij, ial, irr) + varphi_a*(1d0-varphi_r)*phi_X(ij, ix, 0)
                phi_a(ij, iar, irr) = phi_a(ij, iar, irr) + (1d0-varphi_a)*(1d0-varphi_r)*phi_X(ij, ix, 0)
            enddo
        else

            ! iterate over todays cash on hand
            do ix = 0, NX
                do ir = 0, NAR

                    ! interpolate asset decision
                    call linint_Grow(a_plus(ij, ix, ir), a_l, a_u, a_grow, NA, ial, iar, varphi_a)

                    ! restrict values to grid just in case
                    ial = min(ial, NA)
                    iar = min(iar, NA)
                    varphi_a = min(varphi_a, 1d0)

                    ! get end of period asset distribution
                    phi_a(ij, ial, ir) = phi_a(ij, ial, ir) + varphi_a*phi_X(ij, ix, ir)
                    phi_a(ij, iar, ir) = phi_a(ij, iar, ir) + (1d0-varphi_a)*phi_X(ij, ix, ir)
                enddo
            enddo
        endif

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ix, ia, iw, ir
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)
        real*8 :: var_c(JJ), var_y(JJ), var_yp(JJ), var_a(JJ), var_o(JJ)

        ! calculate cohort averages
        c_coh(:)  = 0d0
        y_coh(:)  = 0d0
        yp_coh(:) = 0d0
        a_coh(:)  = 0d0
        o_coh(:)  = 0d0

        do ij = 1, JR-1

            do iw = 1, NW
                y_coh(ij) = y_coh(ij) + (w*eff(ij)*zeta(iw) + pen(ij))*dist_zeta(iw)
            enddo

            yp_coh(ij) = 0d0

            do ix = 0, NX
                c_coh(ij) = c_coh(ij) + c(ij, ix, 0)*phi_X(ij, ix, 0)
            enddo

            if(ij > 1)then
                do ia = 0, NA
                    a_coh(ij) = a_coh(ij) + a(ia)*phi_a(ij-1, ia, 0)
                    o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia, 0)*phi_a(ij-1, ia, 0)
                enddo
            endif
        enddo

        do ij = JR, JJ

            y_coh(ij) = y_coh(ij) + pen(ij)

            do ir = 0, NAR
                do ix = 0, NX
                    yp_coh(ij) = yp_coh(ij)+ ar(ir)/p_a*phi_X(ij, ix, ir)
                    c_coh(ij)  = c_coh(ij) + c(ij, ix, ir)*phi_X(ij, ix, ir)
                enddo

                do ia = 0, NA
                    a_coh(ij) = a_coh(ij) + a(ia)*phi_a(ij-1, ia, ir)
                    o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia, ir)*phi_a(ij-1, ia, ir)
                enddo
            enddo
        enddo

        ! calculate variances
        var_c  = 0d0
        var_y  = 0d0
        var_yp = 0d0
        var_a  = 0d0
        var_o  = 0d0

        do ij = 1, JR-1

            do iw = 1, NW
                var_y(ij) = var_y(ij) + (w*eff(ij)*zeta(iw))**2*dist_zeta(iw)
            enddo

            var_yp(ij) = 0d0

            do ix = 0, NX
                var_c(ij) = var_c(ij) + c(ij, ix, 0)**2*phi_X(ij, ix, 0)
            enddo

            do ia = 0, NA
                if(ij > 1)then
                    var_a(ij) = var_a(ij) + a(ia)**2*phi_a(ij-1, ia, 0)
                    var_o(ij) = var_o(ij) + omega_plus(ij-1, ia, 0)**2*phi_a(ij-1, ia, 0)
                endif
            enddo
        enddo

        do ij = JR, JJ

            var_y(ij) = pen(ij)**2

            do ir = 0, NAR
                do ix = 0, NX
                    var_yp(ij) = var_yp(ij)+ (ar(ir)/p_a)**2*phi_X(ij, ix, ir)
                    var_c(ij)  = var_c(ij) + c(ij, ix, ir)**2*phi_X(ij, ix, ir)
                enddo

                do ia = 0, NA
                    var_a(ij) = var_a(ij) + a(ia)**2*phi_a(ij-1, ia, ir)
                    var_o(ij) = var_o(ij) + omega_plus(ij-1, ia, ir)**2*phi_a(ij-1, ia, ir)
                enddo
            enddo
        enddo

        var_c  = var_c - c_coh**2
        var_y  = var_y - y_coh**2
        var_yp = var_yp- yp_coh**2
        var_a  = var_a - a_coh**2
        var_o  = var_o - o_coh**2

        ! get age dependent variance of eta
        sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

        ! calculate age specific expectations and variance of exp(eta)
        mu_exp = exp(0.5d0*sigma_eta)
        sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

        ! add variance effects
        var_y  = mu_exp**2*var_y + sigma_exp*y_coh**2 + sigma_exp*var_y
        var_c  = mu_exp**2*var_c + sigma_exp*c_coh**2 + sigma_exp*var_c
        var_yp = mu_exp**2*var_yp+ sigma_exp*yp_coh**2+ sigma_exp*var_yp
        var_a  = mu_exp**2*var_a + sigma_exp*a_coh**2 + sigma_exp*var_a

        ! add level effect to averages
        y_coh  = mu_exp*y_coh
        c_coh  = mu_exp*c_coh
        yp_coh = mu_exp*yp_coh
        a_coh  = mu_exp*a_coh

        ! calculate coefficients of variation
        cv_y  = sqrt(max(var_y, 0d0))/max(y_coh, 1d-10)
        cv_c  = sqrt(max(var_c, 0d0))/max(c_coh, 1d-10)
        cv_yp = sqrt(max(var_yp,0d0))/max(yp_coh,1d-10)
        cv_a  = sqrt(max(var_a, 0d0))/max(a_coh, 1d-10)
        cv_o  = sqrt(max(var_o, 0d0))/max(o_coh, 1d-10)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ixmax(JJ), iamax(JJ), irmax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_X(ixmax)
        call check_grid_a(iamax)
        call check_grid_ar(irmax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        write(21, '(a,a,a)')' IJ      CONS    INCOME       Y_P    ASSETS', &
            '     OMEGA     CV(C)     CV(Y)    CV(YP)     CV(A)     CV(O)', &
            '     IXMAX     IAMAX     IRMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,3i10)')ages(ij), c_coh(ij), y_coh(ij), yp_coh(ij), a_coh(ij), &
                    o_coh(ij), cv_c(ij), cv_y(ij), cv_yp(ij), cv_a(ij), cv_o(ij), &
                    ixmax(ij), iamax(ij), irmax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages), y_coh, legend='Labor Income (Mean)')
        call plot(dble(ages), yp_coh, legend='Annuity Payments (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income/Annuities', ylim=(/0d0, 4d0/))

        call plot(dble(ages(2:JJ)), o_coh(2:JJ))
        call execplot(xlabel='Age j', ylabel='Portfolio Share', ylim=(/0d0, 1d0/))

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Liquid Assets')

        call plot(a, omegar_plus)
        call execplot(xlabel='Wealth a', ylabel='Fraction Annuitized omega_r', &
            xlim=(/0d0, 250d0/), ylim=(/0d0, 1d0/))

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, ir

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do ir = 0, NAR
                    if(phi_a(ij, ia, ir) > 1d-8)iamax(ij) = ia
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none
        integer :: ixmax(JJ), ij, ix, ir

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do ir = 0, NAR
                    if(phi_X(ij, ix, ir) > 1d-8)ixmax(ij) = ix
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_ar(irmax)

        implicit none
        integer :: irmax(JJ), ij, ix, ir

        irmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do ir = 0, NAR
                    if(phi_X(ij, ix, ir) > 1d-8)irmax(ij) = ir
                enddo
            enddo
        enddo

    end subroutine

end program
