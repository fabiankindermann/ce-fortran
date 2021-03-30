!##############################################################################
! PROGRAM PortfolioChoiceHousing
!
! ## Portfolio choice with housing
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog10_22m.f90"

program PortfolioChoiceHousing

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

        ! initialize retirement savings grid
        call grid_Cons_Grow(ah, ah_l, ah_u, ah_grow)

        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = min(w*minval(eff(1:JR-1))*minval(eps(:))*zeta(1), pen(JR))
        X_u = (1d0 + (r_f + mu_r + maxval(vtheta(:))))*a_u + &
                    max(w*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW), pen(JR))
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! open files
        open(21, file='output.out')

        ! increase number of iterations for rootfinder
        call setiter_root(1000)

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ix, ia, ih

        ! solve the last period's problem
        if(nu > 0d0)then

            ! the portfolio choice problem
            do ia = 1, NA
                do ih = 0, NAH
                    call solve_portfolio(JJ, ia, ih)
                enddo
            enddo
            omega_plus(JJ, 0, :) = omega_plus(JJ, 1, :)

            call interpolate(JJ)

            ! the retirement account investment problem
            do ia = 1, NA
                call solve_real(JJ, ia)
            enddo
            omegah_plus(JJ, 0) = omegah_plus(JJ, 1)

            call interpolate_RHS(JJ)

            ! the consumption-savings problems
            do ix = 0, NX
                do ih = 0, NAH
                    call solve_consumption(JJ, ix, ih)
                enddo
            enddo

        else
            omega_plus(JJ, :, :) = 0d0
            omegah_plus(JJ, :) = 0d0
            do ix = 0, NX
                do ih = 0, NAH
                    a_plus(JJ, ix, ih) = 0d0
                    c(JJ, ix, ih) = X(ix)
                    V(JJ, ix, ih) = valuefunc(0d0, c(JJ, ix, ih), ah(ih), JJ)
                enddo
            enddo
        endif

        write(*,'(a,i3,a)')'Age: ',JJ,' DONE!'

        ! now household problem during working life
        do ij = JJ-1, 1, -1

            ! the portfolio choice problem
            do ia = 1, NA
                do ih = 0, NAH
                    call solve_portfolio(ij, ia, ih)
                enddo
            enddo
            omega_plus(ij, 0, :) = omega_plus(ij, 1, :)

            ! interpolate individual RHS and temporary value function
            call interpolate(ij)

            ! the retirement account investment problem
            do ia = 1, NA
                call solve_real(ij, ia)
            enddo
            omegah_plus(ij, 0) = omegah_plus(ij, 1)

            ! interpolate individual RHS and temporary value function
            call interpolate_RHS(ij)

            ! the consumption-savings problems
            do ix = 0, NX
                do ih = 0, NAH
                    call solve_consumption(ij, ix, ih)
                enddo
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia, ih)

        implicit none
        integer, intent(in) :: ij, ia, ih
        real*8 :: x_in, port0, port1, tolerance
        real*8 :: x_l, x_r, x_new, f_l, f_new
        integer :: it
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        ih_com = ih

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

        ! use intermediate value theorem
        if(port0*port1 > 0d0)then
            if(abs(port0) > abs(port1))then
                omega_plus(ij, ia, ih) = 1d0
            else
                omega_plus(ij, ia, ih) = 0d0
            endif
            return
        else

            ! get order of magnitude of foc
            tolerance = 1d-5*abs(port0-port1)
            tolerance = min(tolerance, 1d-8)
            call settol_root(tolerance)

            ! get best guess for the root of foc_port
            if(ia > 1)then
                x_in = omega_plus(ij, ia-1, ih)
            else
                x_in = -port0/(port1-port0)
            endif
            check = .false.

            ! solve the household problem using rootfinding
            call fzero(x_in, foc_port, check)

            ! double safety net
            if(check .or. x_in < 0d0 .or. x_in > 1d0)then

                ! perform a bisection search
                check = .false.
                x_l = 0d0
                x_r = 1d0
                f_l = port0
                it = 1
                do
                    ! do the updating step
                    x_new = (x_l + x_r)/2d0
                    f_new = foc_port(x_new)

                    ! check for convergence
                    if(abs(x_new-x_l) < 1d-8)then
                        x_in = x_new
                        exit
                    endif

                    ! check for maximum iterations
                    if(it == 1000)then
                        check = .true.
                        x_in = x_new
                        exit
                    endif

                    ! updating step
                    if(f_l*f_new <= 0d0)then
                        x_r = x_new
                    else
                        x_l = x_new
                        f_l = f_new
                    endif

                    it = it + 1
                enddo

                ! if anything goes wrong then, print error message
                if(check .or. x_in < 0d0 .or. x_in > 1d0)then
                    write(*,'(a, 3i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia, ih
                endif
            endif

            omega_plus(ij, ia, ih) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ix, ih)

        implicit none
        integer, intent(in) :: ij, ix, ih
        real*8 :: x_in
        logical :: check

        ! determine decision for zero cash-on-hand
        if(X(ix) < 1d-10)then
            a_plus(ij, ix, ih) = 0d0
            c(ij, ix, ih) = 0d0
            V(ij, ix, ih) = valuefunc(0d0, 0d0, ah(ih), ij)
            return
        endif

        ! set up communication variables
        ij_com = ij
        ix_com = ix
        ih_com = ih

        ! get best initial guess from future period
        if(ij < JJ)then
            x_in = a_plus(ij+1, ix, ih)
        else
            x_in = X(ix)/2d0
        endif
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 3i4)')'ERROR IN ROOTFINDING CONS : ', ij, ix, ih

        ! check for borrowing constraint
        if(x_in < 0d0)then
            x_in = 0d0
            cons_com = X(ix) - (1d0-xi)*hmin
        endif

        ! copy decisions
        a_plus(ij, ix, ih) = x_in
        c(ij, ix, ih) = cons_com
        V(ij, ix, ih) = valuefunc(x_in, cons_com, ah(ih), ij)

    end subroutine


    ! solve the real estate investment problem
    subroutine solve_real(ij, ia)

        implicit none
        integer, intent(in) :: ij, ia
        integer :: it
        real*8 :: x_in, real0, real1
        real*8 :: x_l, x_r, x_new, f_l, f_new
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia

        ! check for corner solutions
        real0 = foc_real(0d0)
        real1 = foc_real(1d0)

        ! use intermediate value theorem
        if(real0*real1 > 0d0)then
            if(abs(real0) > abs(real1))then
                omegah_plus(ij, ia) = 1d0
            else
                omegah_plus(ij, ia) = 0d0
            endif
            return
        else

            ! get best guess for the root of foc_real
            if(ia > 1)then
                x_in = omegah_plus(ij, ia-1)
            else
                x_in = -real0/(real1-real0)
            endif
            check = .false.

            ! solve the household problem using rootfinding
            call fzero(x_in, foc_real, check)

            ! double safety net
            if(check .or. x_in < 0d0 .or. x_in > 1d0)then

                ! perform a bisection search
                check = .false.
                x_l = 0d0
                x_r = 1d0
                f_l = real0
                it = 1
                do
                    ! do the updating step
                    x_new = (x_l + x_r)/2d0
                    f_new = foc_real(x_new)

                    ! check for convergence
                    if(abs(x_new-x_l) < 1d-8)then
                        x_in = x_new
                        exit
                    endif

                    ! check for maximum iterations
                    if(it == 1000)then
                        check = .true.
                        x_in = x_new
                        exit
                    endif

                    ! updating step
                    if(f_l*f_new <= 0d0)then
                        x_r = x_new
                    else
                        x_l = x_new
                        f_l = f_new
                    endif

                    it = it + 1
                enddo

                ! if anything goes wrong then, print error message
                if(check .or. x_in < 0d0 .or. x_in > 1d0)then
                    write(*,'(a, 2i4)')'ERROR IN ROOTFINDING REAL : ', ij, ia
                endif
            endif

            omegah_plus(ij, ia) = x_in
        endif

    end subroutine


    ! for calculating the values of Q, Q_al, and Q_ah to optimize real estate holdings
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, ih, iw, isr
        real*8 :: X_p, c_p, varphi_X, varphi_h, dist, EV, R_port, beq, ah_p
        integer :: ixl, ixr, ihl, ihr

        Q_al(ij, :, :) = 0d0
        Q_ah(ij, :, :) = 0d0
        Q(ij, :, :) = 0d0

        do ia = 0, NA
            do ih = 0, NAH

                ! when agent is in last year
                if(ij == JJ)then

                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, ih)*(mu_r + vtheta(isr))
                        if(a(ia)-xi*ah(ih) < 0d0)R_port = R_port + load

                        ! get future house size
                        ah_p = ah(ih)

                        ! get bequest
                        beq = R_port*(a(ia)-xi*ah(ih)) + (1d0-delta_h)*ah(ih)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! get derivatives Q_al and Q_ah as well as Q itself
                        Q_al(ij, ia, ih) = Q_al(ij, ia, ih) + dist*R_port*(1d0-psi(ij+1))*margu_b(beq)
                        Q_ah(ij, ia, ih) = Q_ah(ij, ia, ih) + dist* &
                            (1d0-delta_h-xi*R_port)*(1d0-psi(ij+1))*margu_b(beq)
                        Q(ij, ia, ih)    = Q(ij, ia, ih) + dist*(1d0-psi(ij+1))*nu*max(beq, 1d-10)**egam/egam
                    enddo

                ! when agent retired in the next period
                elseif(ij >= JR-1)then

                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, ih)*(mu_r + vtheta(isr))
                        if(a(ia)-xi*ah(ih) < 0d0)R_port = R_port + load
                        X_p = R_port*(a(ia)-xi*ah(ih)) + pen(ij+1) + (1d0-delta_h)*ah(ih)
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get future house size
                        ah_p = ah(ih)
                        call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

                        ! get bequest
                        beq = R_port*(a(ia)-xi*ah(ih)) + (1d0-delta_h)*ah(ih)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! calculate consumption and FOC
                        c_p = varphi_X      *c(ij+1, ixl, ih) + &
                              (1d0-varphi_X)*c(ij+1, ixr, ih)
                        c_p = max(c_p, 1d-10)

                        ! get future value function
                        EV = varphi_X      *(egam*V(ij+1, ixl, ih))**(1d0/egam) + &
                             (1d0-varphi_X)*(egam*V(ij+1, ixr, ih))**(1d0/egam)
                        EV = max(EV, 1d-10)

                        ! get derivatives Q_al and Q_ah as well as Q itself
                        Q_al(ij, ia, ih) = Q_al(ij, ia, ih) + dist*R_port*&
                            (psi(ij+1)*margu_c(c_p) + (1d0-psi(ij+1))*margu_b(beq))
                        Q_ah(ij, ia, ih) = Q_ah(ij, ia, ih) + dist*(psi(ij+1)*margu_h(ah_p) + &
                            (1d0-delta_h-xi*R_port)*(psi(ij+1)*margu_c(c_p) + (1d0-psi(ij+1))*margu_b(beq)))
                        Q(ij, ia, ih)    = Q(ij, ia, ih) + dist*(psi(ij+1)*EV**egam/egam + &
                            (1d0-psi(ij+1))*nu*max(beq, 1d-10)**egam/egam)
                    enddo
                else

                    do isr = 1, NSR
                        do iw = 1, NW

                            R_port = 1d0 + r_f + omega_plus(ij, ia, ih)*(mu_r + vtheta(isr))
                            if(a(ia)-xi*ah(ih) < 0d0)R_port = R_port + load
                            X_p = R_port*(a(ia)-xi*ah(ih))/eps(isr) &
                                  + w*eff(ij+1)*zeta(iw) + (1d0-delta_h)*ah(ih)/eps(isr)
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                            ! get future house size
                            ah_p = ah(ih)/eps(isr)
                            call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

                            ! get bequest
                            beq = R_port*(a(ia)-xi*ah(ih)) + (1d0-delta_h)*ah(ih)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_epsvtheta(isr)

                            ! get future interpolated values
                            if(varphi_X <= varphi_h)then
                                c_p   = varphi_X           *c(ij+1, ixl, ihl) + &
                                        (varphi_h-varphi_X)*c(ij+1, ixr, ihl) + &
                                        (1d0-varphi_h)     *c(ij+1, ixr, ihr)
                                EV    = varphi_X           *(egam*V(ij+1, ixl, ihl))**(1d0/egam) + &
                                        (varphi_h-varphi_X)*(egam*V(ij+1, ixr, ihl))**(1d0/egam) + &
                                        (1d0-varphi_h)     *(egam*V(ij+1, ixr, ihr))**(1d0/egam)
                            else
                                c_p   = varphi_h           *c(ij+1, ixl, ihl) + &
                                        (varphi_X-varphi_h)*c(ij+1, ixl, ihr) + &
                                        (1d0-varphi_X)     *c(ij+1, ixr, ihr)
                                EV    = varphi_h           *(egam*V(ij+1, ixl, ihl))**(1d0/egam) + &
                                        (varphi_X-varphi_h)*(egam*V(ij+1, ixl, ihr))**(1d0/egam) + &
                                        (1d0-varphi_X)     *(egam*V(ij+1, ixr, ihr))**(1d0/egam)
                            endif
                            c_p = max(c_p, 1d-10)

                            ! get derivatives Q_al and Q_ah as well as Q itself
                            Q_al(ij, ia, ih) = Q_al(ij, ia, ih) + dist*R_port*&
                                (psi(ij+1)*margu_c(eps(isr)*c_p) + (1d0-psi(ij+1))*margu_b(beq))
                            Q_ah(ij, ia, ih) = Q_ah(ij, ia, ih) + dist*(psi(ij+1)*margu_h(eps(isr)*ah_p) + &
                                (1d0-delta_h-xi*R_port)*(psi(ij+1)*margu_c(eps(isr)*c_p) +(1d0-psi(ij+1))*margu_b(beq)))
                            Q(ij, ia, ih)    = Q(ij, ia, ih) + dist*(psi(ij+1)*EV**egam/egam + &
                                (1d0-psi(ij+1))*nu*max(beq, 1d-10)**egam/egam)
                        enddo
                    enddo
                endif

                Q_al(ij, ia, ih) = Q_al(ij, ia, ih)**(-gamma)
                Q_ah(ij, ia, ih) = Q_ah(ij, ia, ih)**(-gamma)
                Q(ij, ia, ih)    = (egam*Q(ij, ia, ih))**(1d0/egam)
            enddo
        enddo

    end subroutine


    ! for calculating RHS of FOC for consumption
    subroutine interpolate_RHS(ij)

        implicit none
        integer, intent(in) :: ij
        real*8 :: Qal_p, Qah_p
        real*8 :: al_p, ah_p, varphi_a, varphi_h
        integer :: ia, ih, ial, iar, ihl, ihr

        do ia = 0, NA

            ! determine future liquid wealth and housing wealth
            al_p = (1d0-omegah_plus(ij, ia))*a(ia)
            ah_p = hmin + omegah_plus(ij, ia)*a(ia)/(1d0-xi)

            ! derive interpolation weights
            call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
            call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

            ! get expected value function Q
            if(varphi_a <= varphi_h)then
                Qal_p     = varphi_a           *Q_al(ij, ial, ihl) + &
                            (varphi_h-varphi_a)*Q_al(ij, iar, ihl) + &
                            (1d0-varphi_h)     *Q_al(ij, iar, ihr)
                Qah_p     = varphi_a           *Q_ah(ij, ial, ihl) + &
                            (varphi_h-varphi_a)*Q_ah(ij, iar, ihl) + &
                            (1d0-varphi_h)     *Q_ah(ij, iar, ihr)
                S(ij, ia) = varphi_a           *Q(ij, ial, ihl) + &
                            (varphi_h-varphi_a)*Q(ij, iar, ihl) + &
                            (1d0-varphi_h)     *Q(ij, iar, ihr)
            else
                Qal_p     = varphi_h           *Q_al(ij, ial, ihl) + &
                            (varphi_a-varphi_h)*Q_al(ij, ial, ihr) + &
                            (1d0-varphi_a)     *Q_al(ij, iar, ihr)
                Qah_p     = varphi_h           *Q_ah(ij, ial, ihl) + &
                            (varphi_a-varphi_h)*Q_ah(ij, ial, ihr) + &
                            (1d0-varphi_a)     *Q_ah(ij, iar, ihr)
                S(ij, ia) = varphi_h           *Q(ij, ial, ihl) + &
                            (varphi_a-varphi_h)*Q(ij, ial, ihr) + &
                            (1d0-varphi_a)     *Q(ij, iar, ihr)
            endif
            Qal_p = max(Qal_p, 1d-10)**(-1d0/gamma)
            Qah_p = max(Qah_p, 1d-10)**(-1d0/gamma)

            ! get right hand side of the first order condition
            RHS(ij, ia) = (1d0-omegah_plus(ij, ia))*Qal_p + omegah_plus(ij, ia)*Qah_p/(1d0-xi)
            RHS(ij, ia) = (beta*RHS(ij, ia)/theta)**(-gamma)
        enddo

    end subroutine


    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij

        ! set distributions to zero
        phi_X(:, :, :) = 0d0
        phi_a(:, :)    = 0d0
        phi_h(:, :, :) = 0d0

        do ij = 1, JJ

            ! get distribution on cash-on-hand grid
            call get_distribution_X(ij)

            ! get distribution after consumption-savings decision
            call get_distribution_a(ij)

            ! get distribution after portolio choice decision
            call get_distribution_h(ij)
        enddo

    end subroutine


    ! to calculate distribution on cash-on-hand grid
    subroutine get_distribution_X(ij)

        implicit none
        integer, intent(in) :: ij
        real*8 :: X_p, R_port, varphi_X, dist
        real*8 :: ah_plus, R_j, varphi_h
        integer :: ia, ih, iw, isr, ixl, ixr, ihl, ihr

        if(ij == 1)then

            ! get initial distribution at age 1 of cash-on-hand
            do iw = 1, NW

                ! get initial cash-on-hand
                X_p = w*eff(1)*zeta(iw) + (1d0-delta_h)*eff(1)

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)
                call linint_Grow(eff(1), ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

                ! get distributional weight
                dist = dist_zeta(iw)

                ! initialize the distribution
                phi_X(1, ixl, ihl) = phi_X(1, ixl, ihl) + dist*varphi_X*varphi_h
                phi_X(1, ixr, ihl) = phi_X(1, ixr, ihl) + dist*(1d0-varphi_X)*varphi_h
                phi_X(1, ixl, ihr) = phi_X(1, ixl, ihr) + dist*varphi_X*(1d0-varphi_h)
                phi_X(1, ixr, ihr) = phi_X(1, ixr, ihr) + dist*(1d0-varphi_X)*(1d0-varphi_h)
            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays wealth distribution
            do ia = 0, NA
                do ih = 0, NAH

                    do iw = 1, NW
                        do isr = 1, NSR

                            ! get today's cash-on-hand
                            R_port = 1d0 + (r_f + omega_plus(ij-1, ia, ih)*(mu_r + vtheta(isr)))
                            if(a(ia)-xi*ah(ih) < 0d0)R_port = R_port + load
                            X_p = R_port*(a(ia)-xi*ah(ih))/eps(isr) + w*eff(ij)*zeta(iw) + (1d0-delta_h)*ah(ih)/eps(isr)

                            ! derive interpolation weights
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                            ! get future balance of retirement account and interpolate
                            ah_plus = ah(ih)/eps(isr)
                            call linint_Grow(ah_plus, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_epsvtheta(isr)

                            ! distribute on today's state space
                            phi_X(ij, ixl, ihl) = phi_X(ij, ixl, ihl) + &
                                        dist*varphi_X*varphi_h*phi_h(ij-1, ia, ih)
                            phi_X(ij, ixr, ihl) = phi_X(ij, ixr, ihl) + &
                                        dist*(1d0-varphi_X)*varphi_h*phi_h(ij-1, ia, ih)
                            phi_X(ij, ixl, ihr) = phi_X(ij, ixl, ihr) + &
                                        dist*varphi_X*(1d0-varphi_h)*phi_h(ij-1, ia, ih)
                            phi_X(ij, ixr, ihr) = phi_X(ij, ixr, ihr) + &
                                        dist*(1d0-varphi_X)*(1d0-varphi_h)*phi_h(ij-1, ia, ih)
                        enddo
                    enddo
                enddo
            enddo

        else

            ! iterate over yesterdays wealth distribution
            do ia = 0, NA
                do ih = 0, NAH

                    do isr = 1, NSR

                        ! get today's cash-on-hand
                        R_port = 1d0 + r_f + omega_plus(ij-1, ia, ih)*(mu_r + vtheta(isr))
                        if(a(ia)-xi*ah(ih) < 0d0)R_port = R_port + load
                        X_p = R_port*(a(ia)-xi*ah(ih)) + pen(ij) + (1d0-delta_h)*ah(ih)

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! distribute on today's state space
                        phi_X(ij, ixl, ih) = phi_X(ij, ixl, ih) + dist*varphi_X*phi_h(ij-1, ia, ih)
                        phi_X(ij, ixr, ih) = phi_X(ij, ixr, ih) + dist*(1d0-varphi_X)*phi_h(ij-1, ia, ih)

                    enddo
                enddo
            enddo
        endif

    end subroutine


    ! to calculate distribution after consumption savings decision
    subroutine get_distribution_a(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ix, ih, ial, iar
        real*8 :: varphi_a

        ! iterate over todays cash on hand and retirement wealth
        do ix = 0, NX
            do ih = 0, NAH

                ! interpolate asset decision
                call linint_Grow(a_plus(ij, ix, ih), a_l, a_u, a_grow, NA, ial, iar, varphi_a)

                ! restrict values to grid just in case
                ial = min(ial, NA)
                iar = min(iar, NA)
                varphi_a = min(varphi_a, 1d0)

                ! get distribution over total assets and retirement wealth
                phi_a(ij, ial) = phi_a(ij, ial) + varphi_a*phi_X(ij, ix, ih)
                phi_a(ij, iar) = phi_a(ij, iar) + (1d0-varphi_a)*phi_X(ij, ix, ih)
            enddo
        enddo

    end subroutine


    ! to calculate distribution after real estate decision
    subroutine get_distribution_h(ij)

        implicit none
        integer, intent(in) :: ij
        real*8 :: al_p, ah_p, varphi_a, varphi_h
        integer :: ia, ial, iar, ihl, ihr

        ! iterate over both wealth dimensions
        do ia = 0, NA

            ! determine future liquid wealth and pre-return value of ret. account
            al_p = (1d0-omegah_plus(ij, ia))*a(ia)
            ah_p = hmin + omegah_plus(ij, ia)*a(ia)/(1d0-xi)

            ! derive interpolation weights
            call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
            call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

            ! get distribution over liquid asset and pre-return ret. wealth
            phi_h(ij, ial, ihl) = phi_h(ij, ial, ihl) + varphi_a*varphi_h*phi_a(ij, ia)
            phi_h(ij, iar, ihl) = phi_h(ij, iar, ihl) + (1d0-varphi_a)*varphi_h*phi_a(ij, ia)
            phi_h(ij, ial, ihr) = phi_h(ij, ial, ihr) + varphi_a*(1d0-varphi_h)*phi_a(ij, ia)
            phi_h(ij, iar, ihr) = phi_h(ij, iar, ihr) + (1d0-varphi_a)*(1d0-varphi_h)*phi_a(ij, ia)
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ix, ia, iw, ih
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)
        real*8 :: vah_c(JJ), vah_y(JJ), vah_yp(JJ), vah_a(JJ), vah_ah(JJ), vah_o(JJ)

        ! calculate cohort averages
        c_coh(:)  = 0d0
        y_coh(:)  = 0d0
        a_coh(:)  = 0d0
        ah_coh(:) = 0d0
        o_coh(:)  = 0d0

        ah_coh(1) = eff(1)

        do ij = 1, JJ

            do iw = 1, NW
                y_coh(ij) = y_coh(ij) + (w*eff(ij)*zeta(iw) + pen(ij))*dist_zeta(iw)
            enddo

            do ih = 0, NAH
                do ix = 0, NX
                    c_coh(ij) = c_coh(ij) + c(ij, ix, ih)*phi_X(ij, ix, ih)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        a_coh(ij)  = a_coh(ij) + (a(ia)-xi*ah(ih))*phi_h(ij-1, ia, ih)
                        ah_coh(ij) = ah_coh(ij)+ ah(ih)*phi_h(ij-1, ia, ih)
                        o_coh(ij)  = o_coh(ij) + omega_plus(ij-1, ia, ih)*phi_h(ij-1, ia, ih)
                    enddo
                endif
            enddo
        enddo

        ! calculate variances
        vah_c  = 0d0
        vah_y  = 0d0
        vah_yp = 0d0
        vah_a  = 0d0
        vah_ah = 0d0
        vah_o  = 0d0

        do ij = 1, JJ

            do iw = 1, NW
                vah_y(ij) = vah_y(ij) + (w*eff(ij)*zeta(iw) + pen(ij))**2*dist_zeta(iw)
            enddo

            do ih = 0, NAH
                do ix = 0, NX
                    vah_c(ij) = vah_c(ij) + c(ij, ix, ih)**2*phi_X(ij, ix, ih)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        vah_a(ij)  = vah_a(ij) + (a(ia)-xi*ah(ih))**2*phi_h(ij-1, ia, ih)
                        vah_ah(ij) = vah_ah(ij)+ ah(ih)**2*phi_h(ij-1, ia, ih)
                        vah_o(ij)  = vah_o(ij) + omega_plus(ij-1, ia, ih)**2*phi_h(ij-1, ia, ih)
                    enddo
                else
                    vah_a(ij)  = a_coh(ij)**2
                    vah_ah(ij) = ah_coh(ij)**2
                    vah_o(ij)  = o_coh(ij)**2
                endif
            enddo
        enddo

        vah_c  = vah_c - c_coh**2
        vah_y  = vah_y - y_coh**2
        vah_a  = vah_a - a_coh**2
        vah_ah = vah_ah- ah_coh**2
        vah_o  = vah_o - o_coh**2

        ! get age dependent variance of eta
        sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

        ! calculate age specific expectations and variance of exp(eta)
        mu_exp = exp(0.5d0*sigma_eta)
        sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

        ! add variance effects
        vah_y  = mu_exp**2*vah_y + sigma_exp*y_coh**2 + sigma_exp*vah_y
        vah_c  = mu_exp**2*vah_c + sigma_exp*c_coh**2 + sigma_exp*vah_c
        vah_a  = mu_exp**2*vah_a + sigma_exp*a_coh**2 + sigma_exp*vah_a
        vah_ah = mu_exp**2*vah_ah+ sigma_exp*ah_coh**2+ sigma_exp*vah_ah

        ! add level effect to averages
        y_coh  = mu_exp*y_coh
        c_coh  = mu_exp*c_coh
        a_coh  = mu_exp*a_coh
        ah_coh = mu_exp*ah_coh

        ! calculate coefficients of variation
        cv_y  = sqrt(vah_y)/max(y_coh, 1d-10)
        cv_c  = sqrt(vah_c)/max(c_coh, 1d-10)
        cv_a  = sqrt(vah_a)/max(a_coh, 1d-10)
        cv_ah = sqrt(vah_ah)/max(ah_coh, 1d-10)
        cv_o  = sqrt(vah_o)/max(o_coh, 1d-10)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ixmax(JJ), iamax(JJ), ihmax(JJ), ages(JJ)

        ! check for the maximium grid points used
        call check_grid_X(ixmax)
        call check_grid_a(iamax)
        call check_grid_ah(ihmax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ)/)

        write(21, '(a,a,a)')' IJ      CONS    INCOME       A_H       A_L', &
            '     OMEGA     CV(C)     CV(Y)   CV(A_H)     CV(A)     CV(O)', &
            '     IXMAX     IAMAX     IHMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,3i10)')ages(ij), c_coh(ij), y_coh(ij), ah_coh(ij), a_coh(ij), &
                    o_coh(ij), cv_c(ij), cv_y(ij), cv_ah(ij), cv_a(ij), cv_o(ij), &
                    ixmax(ij), iamax(ij), ihmax(ij)
        enddo
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages), ah_coh, legend='Housing (Mean)')
        call plot(dble(ages), y_coh, legend='Labor Income (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income/Annuities', ylim=(/0d0, 5d0/))

        call plot(dble(ages(2:JJ)), o_coh(2:JJ))
        call execplot(xlabel='Age j', ylabel='Portfolio Share', ylim=(/0d0, 1d0/))

        call plot(dble(ages), a_coh, legend='Liquid Assets')
        call plot(dble(ages), ah_coh, legend='Housing Wealth')
        call execplot(xlabel='Age j', ylabel='Assets', ylim=(/-5d0, 40d0/))

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, ih

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                if(phi_a(ij, ia) > 1d-8)iamax(ij) = ia
                do ih = 0, NAH
                    if(phi_h(ij, ia, ih) > 1d-8)iamax(ij) = ia
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none
        integer :: ixmax(JJ), ij, ix, ih

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do ih = 0, NAH
                    if(phi_X(ij, ix, ih) > 1d-8)ixmax(ij) = ix
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_ah(ihmax)

        implicit none
        integer :: ihmax(JJ), ij, ix, ih

        ihmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do ih = 0, NAH
                    if(phi_X(ij, ix, ih) > 1d-8)ihmax(ij) = ih
                enddo
            enddo
        enddo

    end subroutine

end program
