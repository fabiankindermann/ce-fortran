!##############################################################################
! PROGRAM PortfolioChoice
!
! ## Portfolio choice in the life cycle model (with different aggregations)
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
include "prog10_04m.f90"

program PortfolioChoice

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
        real*8 :: temp(NSR, 2)

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

        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = min(w*minval(eff(1:JR-1))*minval(eps(:))*zeta(1), pen(JR))
        X_u = (1d0 + r_f + mu_r + maxval(vtheta(:)))*a_u + &
                    w*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW)
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! open files
        open(21, file='output.out')

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ix, ia

        ! get decisions in last period of life
        omega_plus(JJ, :) = 0d0
        do ix = 0, NX
            a_plus(JJ, ix) = 0d0
            c(JJ, ix) = X(ix)
            V(JJ, ix) = valuefunc(0d0, c(JJ, ix), JJ)
        enddo

        do ij = JJ-1, 1, -1

            ! determine optimal portfolio choice for all others
            do ia = 1, NA
                call solve_portfolio(ij, ia)
            enddo

            ! set omega for zero savings consistent with next gridpoint
            omega_plus(ij, 0) = omega_plus(ij, 1)

            ! interpolate individual RHS and value function
            call interpolate(ij)

            ! determine consumption-savings solution
            do ix = 0, NX
                call solve_consumption(ij, ix)
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia)

        implicit none
        integer, intent(in) :: ij, ia
        real*8 :: x_in, port0, port1, tolerance
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

        ! use intermediate value theorem
        if(port0*port1 > 0d0)then
            if(abs(port0) > abs(port1))then
                omega_plus(ij, ia) = 1d0
            else
                omega_plus(ij, ia) = 0d0
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
            if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia

            omega_plus(ij, ia) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ix)

        implicit none
        integer, intent(in) :: ij, ix
        real*8 :: x_in
        logical :: check

        ! determine decision for zero cash-on-hand
        if(X(ix) < 1d-10)then
            a_plus(ij, ix) = 0d0
            c(ij, ix) = 0d0
            V(ij, ix) = valuefunc(0d0, 0d0, ij)
            return
        endif

        ! set up communication variables
        ij_com = ij
        ix_com = ix

        ! get best initial guess from future period
        x_in = a_plus(ij+1, ix)
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 2i4)')'ERROR IN ROOTFINDING CONS : ', ij, ix

        ! check for borrowing constraint
        if(x_in < 0d0)then
            x_in = 0d0
            cons_com = X(ix)
        endif

        ! copy decisions
        a_plus(ij, ix) = x_in
        c(ij, ix) = cons_com
        V(ij, ix) = valuefunc(x_in, cons_com, ij)

    end subroutine


    ! for calculating the rhs of the first order condition at age ij
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, iw, isr
        real*8 :: X_p, c_p, varphi, dist, EV, R_port
        integer :: ixl, ixr

        RHS(ij, :) = 0d0
        Q(ij, :) = 0d0

        do ia = 0, NA

            ! case agent is retired tomorrow
            if(ij >= JR-1)then

                do isr = 1, NSR

                    ! get return on the portfolio
                    R_port = 1d0 + r_f + omega_plus(ij, ia)*(mu_r + vtheta(isr))

                    ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                    X_p = R_port*a(ia) + pen(ij+1)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_epsvtheta(isr)

                    ! get future consumption value
                    c_p = max(varphi*c(ij+1, ixl) + (1d0-varphi)*c(ij+1, ixr), 1d-10)

                    ! get tomorrow's value function
                    EV = varphi      *(egam*V(ij+1, ixl))**(1d0/egam) + &
                         (1d0-varphi)*(egam*V(ij+1, ixr))**(1d0/egam)

                    ! get RHS of foc and Q
                    RHS(ij, ia) = RHS(ij, ia) + dist*R_port*margu(c_p)
                    Q(ij, ia)   = Q(ij, ia) + dist*EV**egam/egam
                enddo

            ! agent is working
            else
                do iw = 1, NW
                    do isr = 1, NSR

                        ! get return on the portfolio
                        R_port = 1d0 + r_f + omega_plus(ij, ia)*(mu_r + vtheta(isr))

                        ! get tomorrow's cash on hand
                        X_p = R_port*a(ia)/eps(isr) + w*eff(ij+1)*zeta(iw)

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                        ! get distributional weight
                        dist = dist_zeta(iw)*dist_epsvtheta(isr)

                        ! get future consumption value
                        c_p = max(varphi*c(ij+1, ixl) + (1d0-varphi)*c(ij+1, ixr), 1d-10)

                        ! get tomorrow's value function
                        EV = varphi      *(egam*V(ij+1, ixl))**(1d0/egam) + &
                             (1d0-varphi)*(egam*V(ij+1, ixr))**(1d0/egam)

                        ! get RHS of foc and Q
                        RHS(ij, ia) = RHS(ij, ia) + dist*R_port*margu(eps(isr)*c_p)
                        Q(ij, ia)   = Q(ij, ia) + dist*(eps(isr)*EV)**egam/egam
                    enddo
                enddo
            endif

            RHS(ij, ia) = (beta*psi(ij+1)*RHS(ij, ia))**(-gamma)
            Q(ij, ia)   = (egam*Q(ij, ia))**(1d0/egam)
        enddo

    end subroutine


    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij

        ! set distributions to zero
        phi_X(:, :) = 0d0
        phi_a(:, :) = 0d0

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
        integer, intent(in) :: ij
        integer :: ia, iw, isr, ixl, ixr
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
                phi_X(1, ixl) = phi_X(1, ixl) + dist*varphi
                phi_X(1, ixr) = phi_X(1, ixr) + dist*(1d0-varphi)
            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays asset distribution
            do ia = 0, NA

                ! iterate over current shocks
                do iw = 1, NW
                    do isr = 1, NSR

                        ! get today's cash-on-hand
                        R_port = 1d0 + r_f + omega_plus(ij-1, ia)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia)/eps(isr) + w*eff(ij)*zeta(iw)

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                        ! get distributional weight
                        dist = dist_zeta(iw)*dist_epsvtheta(isr)

                        phi_X(ij, ixl) = phi_X(ij, ixl) + dist*varphi*phi_a(ij-1, ia)
                        phi_X(ij, ixr) = phi_X(ij, ixr) + dist*(1d0-varphi)*phi_a(ij-1, ia)
                    enddo
                enddo
            enddo

        else

            ! iterate over yesterdays asset distribution
            do ia = 0, NA

                ! iterate over current shocks
                do isr = 1, NSR

                    ! get today's cash-on-hand
                    R_port = 1d0 + r_f + omega_plus(ij-1, ia)*(mu_r + vtheta(isr))
                    X_p = R_port*a(ia) + pen(ij)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_epsvtheta(isr)

                    phi_X(ij, ixl) = phi_X(ij, ixl) + dist*varphi*phi_a(ij-1, ia)
                    phi_X(ij, ixr) = phi_X(ij, ixr) + dist*(1d0-varphi)*phi_a(ij-1, ia)
                enddo
            enddo
        endif

    end subroutine


    ! to calculate the end of period asset distribution
    subroutine get_distribution_a(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ix, ial, iar
        real*8 :: varphi

        ! iterate over todays cash on hand
        do ix = 0, NX

            ! interpolate asset decision
            call linint_Grow(a_plus(ij, ix), a_l, a_u, a_grow, NA, ial, iar, varphi)

            ! restrict values to grid just in case
            ial = min(ial, NA)
            iar = min(iar, NA)
            varphi = min(varphi, 1d0)

            ! get end of period asset distribution
            phi_a(ij, ial) = phi_a(ij, ial) + varphi*phi_X(ij, ix)
            phi_a(ij, iar) = phi_a(ij, iar) + (1d0-varphi)*phi_X(ij, ix)
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ix, ia, ie, iw
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)
        real*8 :: var_c(JJ), var_a(JJ), var_y(JJ), var_o(JJ)

        ! generate eta distribution if not analytical calculation
        if(.not. analytical)call generate_eta()

        ! calculate cohort averages
        c_coh(:) = 0d0
        y_coh(:) = 0d0
        a_coh(:) = 0d0
        o_coh(:) = 0d0

        ! analytical approach or not
        if(analytical)then

            do ij = 1, JJ

                do iw = 1, NW
                    y_coh(ij) = y_coh(ij) + (w*eff(ij)*zeta(iw) + pen(ij))*dist_zeta(iw)
                enddo

                do ix = 0, NX
                    c_coh(ij) = c_coh(ij) + c(ij, ix)*phi_X(ij, ix)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        a_coh(ij) = a_coh(ij) + a(ia)*phi_a(ij-1, ia)
                        o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia)*phi_a(ij-1, ia)
                    enddo
                endif
            enddo

        else

            do ij = 1, JJ
                do ie = 0, NE

                    do iw = 1, NW
                        y_coh(ij) = y_coh(ij) + (w*eff(ij)*zeta(iw) + pen(ij))*eta(ij, ie) &
                                             *dist_zeta(iw)*phi_e(ij, ie)
                    enddo

                    do ix = 0, NX
                        c_coh(ij) = c_coh(ij) + c(ij, ix)*eta(ij, ie)*phi_X(ij, ix)*phi_e(ij, ie)
                    enddo

                    if(ij > 1)then
                        do ia = 0, NA
                            a_coh(ij) = a_coh(ij) + a(ia)*eta(ij, ie)*phi_a(ij-1, ia)*phi_e(ij, ie)
                            o_coh(ij) = o_coh(ij) + omega_plus(ij-1, ia)*phi_a(ij-1, ia)*phi_e(ij, ie)
                        enddo
                    endif
                enddo
            enddo

        endif

        ! calculate variances
        var_c = 0d0
        var_y = 0d0
        var_a = 0d0
        var_o = 0d0

        ! analytical approach or not
        if(analytical)then

            do ij = 1, JJ

                do iw = 1, NW
                    var_y(ij) = var_y(ij) + ((w*eff(ij)*zeta(iw) + pen(ij)))**2*dist_zeta(iw)
                enddo

                do ix = 0, NX
                    var_c(ij) = var_c(ij) + c(ij, ix)**2*phi_X(ij, ix)
                enddo

                do ia = 0, NA
                    if(ij > 1)then
                        var_a(ij) = var_a(ij) + a(ia)**2*phi_a(ij-1, ia)
                        var_o(ij) = var_o(ij) + omega_plus(ij-1, ia)**2*phi_a(ij-1, ia)
                    endif
                enddo
            enddo

        else

            do ij = 1, JJ
                do ie = 0, NE

                    do iw = 1, NW
                        var_y(ij) = var_y(ij) + ((w*eff(ij)*zeta(iw) + pen(ij))*eta(ij, ie))**2 &
                                             *dist_zeta(iw)*phi_e(ij, ie)
                    enddo

                    do ix = 0, NX
                        var_c(ij) = var_c(ij) + (c(ij, ix)*eta(ij, ie))**2*phi_X(ij, ix)*phi_e(ij, ie)
                    enddo

                    if(ij > 1)then
                        do ia = 0, NA
                            var_a(ij) = var_a(ij) + (a(ia)*eta(ij, ie))**2*phi_a(ij-1, ia)*phi_e(ij, ie)
                            var_o(ij) = var_o(ij) + omega_plus(ij-1, ia)**2*phi_a(ij-1, ia)*phi_e(ij, ie)
                        enddo
                    endif
                enddo
            enddo
        endif

        var_c = var_c - c_coh**2
        var_y = var_y - y_coh**2
        var_a = var_a - a_coh**2
        var_o = var_o - o_coh**2

        ! in case of analytical, add extra eta effect
        if(analytical)then

            ! get age dependent variance of eta
            sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

            ! calculate age specific expectations and variance of exp(eta)
            mu_exp = exp(0.5d0*sigma_eta)
            sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

            ! add variance effects
            var_y = mu_exp**2*var_y + sigma_exp*y_coh**2 + sigma_exp*var_y
            var_c = mu_exp**2*var_c + sigma_exp*c_coh**2 + sigma_exp*var_c
            var_a = mu_exp**2*var_a + sigma_exp*a_coh**2 + sigma_exp*var_a

            ! add level effect to averages
            y_coh = mu_exp*y_coh
            c_coh = mu_exp*c_coh
            a_coh = mu_exp*a_coh
        endif

        ! calculate coefficients of variation
        cv_y = sqrt(max(var_y, 0d0))/max(y_coh, 1d-10)
        cv_c = sqrt(max(var_c, 0d0))/max(c_coh, 1d-10)
        cv_a = sqrt(max(var_a, 0d0))/max(a_coh, 1d-10)
        cv_o = sqrt(max(var_o, 0d0))/max(o_coh, 1d-10)

        if(.not. analytical)call calculate_quantiles

    end subroutine


    ! generates the eta distribution
    subroutine generate_eta()

        implicit none
        integer :: ij, iel, ier, ie, isr
        real*8 :: varphi, eta_temp

        ! set bounds and grid for working ages
        eta_l(1)  = 0d0
        eta_u(1)  = 0d0
        eta(1, :) = 0d0

        do ij = 2, JR-1
            eta_l(ij)  = (ij-1)*(minval(log(eps)))
            eta_u(ij)  = (ij-1)*(maxval(log(eps)))
            call grid_Cons_Equi(eta(ij, :), eta_l(ij), eta_u(ij))
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            eta_l(ij)  = eta_l(JR-1)
            eta_u(ij)  = eta_u(JR-1)
            eta(ij, :) = eta(JR-1, :)
        enddo

        phi_e = 0d0

        ! initial distribution at age 1
        phi_e(1, :) = 1d0/dble(NE+1)

        ! iterate over working different years
        do ij = 2, JR-1

            ! last period's etas
            do ie = 0, NE

                ! new innovations
                do isr = 1, NSR

                    ! distribute on new grid
                    eta_temp = eta(ij-1, ie) + log(eps(isr))
                    call linint_Equi(eta_temp, eta_l(ij), eta_u(ij), NE, iel, ier, varphi)
                    phi_e(ij, iel) = phi_e(ij, iel) + dist_epsvtheta(isr)*varphi*phi_e(ij-1, ie)
                    phi_e(ij, ier) = phi_e(ij, ier) + dist_epsvtheta(isr)*(1d0-varphi)*phi_e(ij-1, ie)
                enddo
            enddo
        enddo

        ! no innovations throughout retirement
        do ij = JR, JJ
            phi_e(ij, :) = phi_e(JR-1, :)
        enddo

        ! take exponentials
        eta = exp(eta)

    end subroutine


    ! subroutine to calculate age specific quantiles of the distribution
    subroutine calculate_quantiles()

        use toolbox

        implicit none
        integer :: ij, ie, ia, ic, it, NC
        real*8 :: a_sort((NA+1)*(NE+1)), a_dist((NA+1)*(NE+1)), a_cdist((NA+1)*(NE+1))
        integer :: iorder((NA+1)*(NE+1))
        real*8 :: thresholds(5), quantiles(size(thresholds, 1), JJ), slope, ages(JJ)

        ! define quantile thresholds
        thresholds = (/0.05d0, 0.25d0, 0.50d0, 0.75d0, 0.95d0/)
        quantiles = 0d0

        ! iterate over ages
        do ij = 2, JJ

            a_sort = 0d0
            a_dist = 0d0

            ! copy savings into one-dimensional array
            ic = 1
            do ie = 0, NE
                do ia = 0, NA
                    if(phi_a(ij-1, ia)*phi_e(ij, ie) > 1d-12)then
                        a_sort(ic) = a(ia)*eta(ij, ie)
                        a_dist(ic) = phi_a(ij-1, ia)*phi_e(ij, ie)
                        ic = ic + 1
                    endif
                enddo
            enddo
            NC = ic -1

            ! sort array and distribution
            call sort(a_sort(1:NC), iorder(1:NC))

            ! calculate cumulative distribution (attention ordering)
            a_cdist(1) = a_dist(iorder(1))
            do ic = 2, NC
                a_cdist(ic) = a_cdist(ic-1) + a_dist(iorder(ic))
            enddo

            ! get quantiles
            do it = 1, size(thresholds, 1)
                if(thresholds(it) <= a_cdist(1))then
                    quantiles(it, ij) = a_sort(1)
                else
                    do ic = 2, NC
                        if(thresholds(it) < a_cdist(ic))then
                            slope = (a_sort(ic)-a_sort(ic-1))/(a_cdist(ic)-a_cdist(ic-1))
                            quantiles(it, ij) = a_sort(ic-1) + slope*(thresholds(it)-a_cdist(ic-1))
                            exit
                        elseif(ic == NC)then
                            quantiles(it, ij) = a_sort(NC)
                        endif
                    enddo
                endif
            enddo

        enddo

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        call plot(ages, quantiles(1, :), legend='5%')
        call plot(ages, quantiles(2, :), legend='25%')
        call plot(ages, quantiles(3, :), legend='50%')
        call plot(ages, quantiles(4, :), legend='75%')
        call plot(ages, quantiles(5, :), legend='95%')

        call execplot(title='Quantiles of Asset Distributions', &
            xlabel='Age j', ylabel='Assets')

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ixmax(JJ), iamax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_X(ixmax)
        call check_grid_a(iamax)

        ! set up age variable
        ages = 20 + (/(ij, ij=1,JJ)/)

        write(21, '(a,a)')' IJ      CONS    INCOME    ASSETS', &
            '     OMEGA     CV(C)     CV(Y)     CV(A)     CV(O)     IXMAX     IAMAX'
        do ij = 1, JJ
            write(21,'(i3,8f10.3,2i10)')ages(ij), c_coh(ij), y_coh(ij), a_coh(ij), o_coh(ij), &
                    cv_c(ij), cv_y(ij), cv_a(ij), cv_o(ij), ixmax(ij), iamax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages), y_coh, legend='Labor Income (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income', ylim=(/0d0, 4d0/))

        call plot(dble(ages(2:JJ)), o_coh(2:JJ))
        call execplot(xlabel='Age j', ylabel='Portfolio Share', ylim=(/0d0, 1d0/))

        call plot(dble(ages), a_coh)
        call execplot(xlabel='Age j', ylabel='Assets')

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                if(phi_a(ij, ia) > 1d-8)iamax(ij) = ia
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none
        integer :: ixmax(JJ), ij, ix

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                if(phi_X(ij, ix) > 1d-8)ixmax(ij) = ix
            enddo
        enddo

    end subroutine

end program
