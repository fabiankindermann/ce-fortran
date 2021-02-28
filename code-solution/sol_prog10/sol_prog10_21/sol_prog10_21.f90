!##############################################################################
! PROGRAM PortfolioChoiceAnnuities
!
! ## Portfolio choice with the continous option to annuitize
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
include "sol_prog10_21m.f90"

program PortfolioChoiceRetAccounts

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
        call grid_Cons_Grow(ya, ya_l, ya_u, ya_grow)

        ! calculate the price of an annuity
        p_a(JJ) = 0d0
        do ij = JJ-1, 1, -1
            p_a(ij) = psi(ij+1)/(1d0+r_f)*(1d0+p_a(ij+1))
        enddo
        p_a(:) = (1d0+xi)*p_a(:)

        ! endogenous lower and upper bound of cash-on-hand grid
        X_l = min(w*minval(eff(1:JR-1))*minval(eps(:))*zeta(1), pen(JR))
        X_u = (1d0 + r_f + mu_r + maxval(vtheta(:)))*a_u + ya_u + &
                    w*maxval(eff(1:JR-1))*maxval(eps(:))*zeta(NW)
        call grid_Cons_Grow(X, X_l, X_u, X_grow)

        ! open files
        open(21, file='output.out')

        ! increase number of iterations for rootfinder
        call setiter_root(1500)

    end subroutine


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ij, ix, ia, iy

        ! solve the last period's problem
        if(nu > 0d0)then

            ! the portfolio choice problem
            do ia = 1, NA
                do iy = 0, NY
                    call solve_portfolio(JJ, ia, iy)
                enddo
            enddo
            omega_plus(JJ, 0, :) = omega_plus(JJ, 1, :)

            call interpolate(JJ)

            omegay_plus(JJ, :, :) = 0d0

            call interpolate_RHS(JJ)

            ! the consumption-savings problems
            do ix = 0, NX
                do iy = 0, NY
                    call solve_consumption(JJ, ix, iy)
                enddo
            enddo

        else
            omega_plus(JJ, :, :) = 0d0
            omegay_plus(JJ, :, :) = 0d0
            V_ya(JJ, :, :) = 0d0
            do ix = 0, NX
                a_plus(JJ, ix, :) = 0d0
                c(JJ, ix, :) = X(ix)
                V(JJ, ix, :) = valuefunc(0d0, c(JJ, ix, 0), JJ, 0)
            enddo
        endif

        write(*,'(a,i3,a)')'Age: ',JJ,' DONE!'


        ! solve household problem
        do ij = JJ-1, 1, -1

            ! the portfolio choice problem
            do ia = 1, NA
                do iy = 0, NY
                    call solve_portfolio(ij, ia, iy)
                enddo
            enddo
            omega_plus(ij, 0, :) = omega_plus(ij, 1, :)

            call interpolate(ij)

            ! the retirement account investment problem
            if(ij < JJ)then
                do ia = 1, NA
                    do iy = 0, NY
                        call solve_annuities(ij, ia, iy)
                    enddo
                enddo
                omegay_plus(ij, 0, :) = omegay_plus(ij, 1, :)
            else
                omegay_plus(ij, :, :) = 0d0
            endif

            call interpolate_RHS(ij)

            ! the consumption-savings problems
            do ix = 0, NX
                do iy = 0, NY
                    call solve_consumption(ij, ix, iy)
                enddo
            enddo

            write(*,'(a,i3,a)')'Age: ',ij,' DONE!'
        enddo

    end subroutine


    ! solve the household's portfolio decision
    subroutine solve_portfolio(ij, ia, iy)

        implicit none
        integer, intent(in) :: ij, ia, iy
        real*8 :: x_in, port0, port1, tolerance
        real*8 :: x_l, x_r, x_new, f_l, f_new
        integer :: it
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        ir_com = iy

        ! check for corner solutions
        port0 = foc_port(0d0)
        port1 = foc_port(1d0)

        ! use intermediate value theorem
        if(port0*port1 > 0d0)then
            if(abs(port0) > abs(port1))then
                omega_plus(ij, ia, iy) = 1d0
            else
                omega_plus(ij, ia, iy) = 0d0
            endif
            return
        else

            ! get order of magnitude of foc
            tolerance = 1d-5*abs(port0-port1)
            tolerance = min(tolerance, 1d-8)
            call settol_root(tolerance)

            ! get best guess for the root of foc_port
            if(ia > 1)then
                x_in = omega_plus(ij, ia-1, iy)
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
                    write(*,'(a, 3i4)')'ERROR IN ROOTFINDING PORT : ', ij, ia, iy
                endif
            endif

            omega_plus(ij, ia, iy) = x_in

            ! reset tolerance level to original value
            call settol_root(1d-8)
        endif

    end subroutine


    ! solve the household's consumption-savings decision
    subroutine solve_consumption(ij, ix, iy)

        implicit none
        integer, intent(in) :: ij, ix, iy
        real*8 :: x_in
        logical :: check

        ! determine decision for zero cash-on-hand
        if(X(ix) < 1d-10)then
            a_plus(ij, ix, iy) = 0d0
            c(ij, ix, iy) = 0d0
            V(ij, ix, iy) = valuefunc(0d0, 0d0, ij, iy)
            V_ya(ij, ix, iy) = dV_dya(0d0, ij, iy)
            return
        endif

        ! set up communication variables
        ij_com = ij
        ix_com = ix
        ir_com = iy

        ! get best initial guess from future period
        if(ij < JJ)then
            x_in = a_plus(ij+1, ix, iy)
        else
            x_in = X(ix)/2d0
        endif
        check = .false.

        ! solve the household problem using rootfinding
        call fzero(x_in, foc_cons, check)

        ! write screen output in case of a problem
        if(check)write(*,'(a, 3i4)')'ERROR IN ROOTFINDING CONS : ', ij, ix, iy

        ! check for borrowing constraint
        if(x_in < 0d0)then
            x_in = 0d0
            cons_com = X(ix)
        endif

        ! copy decisions
        a_plus(ij, ix, iy) = x_in
        c(ij, ix, iy) = cons_com
        V(ij, ix, iy) = valuefunc(x_in, cons_com, ij, iy)
        if(ij == JJ)then
            V_ya(ij, ix, iy) = 0d0
        else
            V_ya(ij, ix, iy) = dV_dya(x_in, ij, iy)
        endif

    end subroutine


    ! solve the retirement account accumulation decision
    subroutine solve_annuities(ij, ia, iy)

        implicit none
        integer, intent(in) :: ij, ia, iy
        integer :: it
        real*8 :: x_in, annu0, annu1
        real*8 :: x_l, x_r, x_new, f_l, f_new
        logical :: check

        ! set up communication variables
        ij_com = ij
        ia_com = ia
        ir_com = iy

        ! check for corner solutions
        annu0 = foc_annu(0d0)
        annu1 = foc_annu(1d0)

        ! use intermediate value theorem
        if(annu0*annu1 > 0d0)then
            if(abs(annu0) > abs(annu1))then
                omegay_plus(ij, ia, iy) = 1d0
            else
                omegay_plus(ij, ia, iy) = 0d0
            endif
            return
        else

            ! get best guess for the root of foc_annu
            if(ia > 1)then
                x_in = omegay_plus(ij, ia-1, iy)
            else
                x_in = -annu0/(annu1-annu0)
            endif
            check = .false.

            ! solve the household problem using rootfinding
            call fzero(x_in, foc_annu, check)

            ! double safety net
            if(check .or. x_in < 0d0 .or. x_in > 1d0)then

                ! perform a bisection search
                check = .false.
                x_l = 0d0
                x_r = 1d0
                f_l = annu0
                it = 1
                do
                    ! do the updating step
                    x_new = (x_l + x_r)/2d0
                    f_new = foc_annu(x_new)

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
                    write(*,'(a, 3i4)')'ERROR IN ROOTFINDING ANNU : ', ij, ia, iy
                endif
            endif

            omegay_plus(ij, ia, iy) = x_in
        endif

    end subroutine


    ! for calculating the values of Q, Q_al, and Q_ya to optimize annuity purchases
    subroutine interpolate(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: iw, isr, ia, iy
        real*8 :: X_p, c_p, varphi_X, dist, EV, R_port
        real*8 :: varphi_y, Vya_p, ya_plus
        integer :: ixl, ixr, iyl, iyr

        Q_al(ij, :, :) = 0d0
        Q_ya(ij, :, :) = 0d0
        Q(ij, :, :) = 0d0

        do ia = 0, NA
            do iy = 0, NY

                ! when agent is in last year
                if(ij == JJ)then
                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, iy)*(mu_r + vtheta(isr))

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! get derivatives Q_al and Q_ya as well as Q itself
                        Q_al(ij, ia, iy) = Q_al(ij, ia, iy) + dist*(1d0-psi(ij+1))*R_port*nu*margu(R_port*a(ia))
                        Q(ij, ia, iy)    = Q(ij, ia, iy) + dist*(1d0-psi(ij+1))*nu*max(R_port*a(ia), 1d-10)**egam/egam
                    enddo

                ! when agent retired in the next period
                elseif(ij >= JR-1)then
                    do isr = 1, NSR

                        ! get future cash-on-hand and interpolate
                        R_port = 1d0 + r_f + omega_plus(ij, ia, iy)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia) + pen(ij+1) + ya(iy)
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! get future consumption value
                        c_p = varphi_X      *c(ij+1, ixl, iy) + &
                              (1d0-varphi_X)*c(ij+1, ixr, iy)
                        c_p = max(c_p, 1d-10)

                        ! get derivative of value function
                        if(ij < JJ-1)then
                            Vya_p = varphi_X      *V_ya(ij+1, ixl, iy) + &
                                    (1d0-varphi_X)*V_ya(ij+1, ixr, iy)
                            Vya_p = max(Vya_p, 1d-10)**(-1d0/gamma)
                        else
                            Vya_p = 0d0
                        endif

                        ! get tomorrow's value function
                        EV = varphi_X      *(egam*V(ij+1, ixl, iy))**(1d0/egam) + &
                             (1d0-varphi_X)*(egam*V(ij+1, ixr, iy))**(1d0/egam)

                        ! get derivatives Q_al and Q_ya as well as Q itself
                        Q_al(ij, ia, iy) = Q_al(ij, ia, iy) + dist*R_port*(psi(ij+1)*margu(c_p) + &
                            nu*(1d0-psi(ij+1))*margu(R_port*a(ia)))
                        Q_ya(ij, ia, iy) = Q_ya(ij, ia, iy) + dist*psi(ij+1)*(margu(c_p) + Vya_p)
                        Q(ij, ia, iy)    = Q(ij, ia, iy) + dist*(psi(ij+1)*EV**egam/egam + &
                            (1d0-psi(ij+1))*nu*max(R_port*a(ia), 1d-10)**egam/egam)
                    enddo
                else

                    do isr = 1, NSR
                        do iw = 1, NW

                            ! get future cash-on-hand and interpolate
                            R_port = 1d0 + r_f + omega_plus(ij, ia, iy)*(mu_r + vtheta(isr))
                            X_p = R_port*a(ia)/eps(isr) + w*eff(ij+1)*zeta(iw) + ya(iy)
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                            ! get future balance of annuity value and interpolate
                            ya_plus = ya(iy)/eps(isr)
                            call linint_Grow(ya_plus, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_epsvtheta(isr)

                            ! get future interpolated values
                            if(varphi_X <= varphi_y)then
                                c_p   = varphi_X           *c(ij+1, ixl, iyl) + &
                                        (varphi_y-varphi_X)*c(ij+1, ixr, iyl) + &
                                        (1d0-varphi_y)     *c(ij+1, ixr, iyr)
                                Vya_p = varphi_X           *V_ya(ij+1, ixl, iyl) + &
                                        (varphi_y-varphi_X)*V_ya(ij+1, ixr, iyl) + &
                                        (1d0-varphi_y)     *V_ya(ij+1, ixr, iyr)
                                EV    = varphi_X           *(egam*V(ij+1, ixl, iyl))**(1d0/egam) + &
                                        (varphi_y-varphi_X)*(egam*V(ij+1, ixr, iyl))**(1d0/egam) + &
                                        (1d0-varphi_y)     *(egam*V(ij+1, ixr, iyr))**(1d0/egam)
                            else
                                c_p   = varphi_y           *c(ij+1, ixl, iyl) + &
                                        (varphi_X-varphi_y)*c(ij+1, ixl, iyr) + &
                                        (1d0-varphi_X)     *c(ij+1, ixr, iyr)
                                Vya_p = varphi_y           *V_ya(ij+1, ixl, iyl) + &
                                        (varphi_X-varphi_y)*V_ya(ij+1, ixl, iyr) + &
                                        (1d0-varphi_X)     *V_ya(ij+1, ixr, iyr)
                                EV    = varphi_y           *(egam*V(ij+1, ixl, iyl))**(1d0/egam) + &
                                        (varphi_X-varphi_y)*(egam*V(ij+1, ixl, iyr))**(1d0/egam) + &
                                        (1d0-varphi_X)     *(egam*V(ij+1, ixr, iyr))**(1d0/egam)
                            endif
                            c_p = max(c_p, 1d-10)
                            Vya_p = max(Vya_p, 1d-10)**(-1d0/gamma)

                            ! get derivatives Q_al and Q_ya as well as Q itself
                            Q_al(ij, ia, iy) = Q_al(ij, ia, iy) + dist*R_port*(psi(ij+1)*margu(eps(isr)*c_p) + &
                                nu*(1d0-psi(ij+1))*margu(R_port*a(ia)))
                            Q_ya(ij, ia, iy) = Q_ya(ij, ia, iy) + dist*psi(ij+1)*(margu(eps(isr)*c_p) &
                                + eps(isr)**(-1d0/gamma)*Vya_p)
                            Q(ij, ia, iy)    = Q(ij, ia, iy) + dist*(psi(ij+1)*(eps(isr)*EV)**egam/egam + &
                                nu*(1d0-psi(ij+1))*max(R_port*a(ia), 1d-10)**egam/egam)
                        enddo
                    enddo
                endif

                Q_al(ij, ia, iy) = Q_al(ij, ia, iy)**(-gamma)
                if(ij < JJ)Q_ya(ij, ia, iy) = Q_ya(ij, ia, iy)**(-gamma)
                Q(ij, ia, iy)    = (egam*Q(ij, ia, iy))**(1d0/egam)
            enddo
        enddo

    end subroutine


    ! for calculating RHS of FOC for consumption
    subroutine interpolate_RHS(ij)

        implicit none
        integer, intent(in) :: ij
        real*8 :: Qal_p, Qya_p
        real*8 :: al_p, ya_p, varphi_a, varphi_y
        integer :: ia, iy, ial, iar, iyl, iyr

        do ia = 0, NA
            do iy = 0, NY

                ! determine future liquid wealth and annuity value
                al_p = (1d0-omegay_plus(ij, ia, iy))*a(ia)
                ya_p = ya(iy) + omegay_plus(ij, ia, iy)*a(ia)/max(p_a(ij), 1d-10)

                ! derive interpolation weights
                call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                call linint_Grow(ya_p, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

                ! get expected value function Q
                if(varphi_a <= varphi_y)then
                    Qal_p         = varphi_a           *Q_al(ij, ial, iyl) + &
                                    (varphi_y-varphi_a)*Q_al(ij, iar, iyl) + &
                                    (1d0-varphi_y)     *Q_al(ij, iar, iyr)
                    Qya_p         = varphi_a           *Q_ya(ij, ial, iyl) + &
                                    (varphi_y-varphi_a)*Q_ya(ij, iar, iyl) + &
                                    (1d0-varphi_y)     *Q_ya(ij, iar, iyr)
                    S(ij, ia, iy) = varphi_a           *Q(ij, ial, iyl) + &
                                    (varphi_y-varphi_a)*Q(ij, iar, iyl) + &
                                    (1d0-varphi_y)     *Q(ij, iar, iyr)
                else
                    Qal_p         = varphi_y           *Q_al(ij, ial, iyl) + &
                                    (varphi_a-varphi_y)*Q_al(ij, ial, iyr) + &
                                    (1d0-varphi_a)     *Q_al(ij, iar, iyr)
                    Qya_p         = varphi_y           *Q_ya(ij, ial, iyl) + &
                                    (varphi_a-varphi_y)*Q_ya(ij, ial, iyr) + &
                                    (1d0-varphi_a)     *Q_ya(ij, iar, iyr)
                    S(ij, ia, iy) = varphi_y           *Q(ij, ial, iyl) + &
                                    (varphi_a-varphi_y)*Q(ij, ial, iyr) + &
                                    (1d0-varphi_a)     *Q(ij, iar, iyr)
                endif
                Qal_p = max(Qal_p, 1d-10)**(-1d0/gamma)
                Qya_p = max(Qya_p, 1d-10)**(-1d0/gamma)

                ! get right hand side of the first order condition
                RHS(ij, ia, iy) = (1d0-omegay_plus(ij, ia, iy))*Qal_p &
                                   + omegay_plus(ij, ia, iy)*Qya_p/max(p_a(ij), 1d-10)
                RHS(ij, ia, iy) = (beta*RHS(ij, ia, iy))**(-gamma)

                ! store the value of Qya_p in Var_temp
                S_ya(ij, ia, iy) = Qya_p**(-gamma)

            enddo
        enddo

    end subroutine


    ! determines the invariant distribution over state space
    subroutine get_distribution()

        implicit none
        integer :: ij

        ! set distributions to zero
        phi_X(:, :, :) = 0d0
        phi_a(:, :, :) = 0d0
        phi_r(:, :, :) = 0d0

        do ij = 1, JJ

            ! get distribution on cash-on-hand grid
            call get_distribution_X(ij)

            ! get distribution after consumption-savings decision
            call get_distribution_a(ij)

            ! get distribution after portolio choice decision
            call get_distribution_ya(ij)
        enddo

    end subroutine


    ! to calculate distribution on cash-on-hand grid
    subroutine get_distribution_X(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ia, iy, iw, isr, ixl, ixr, iyl, iyr
        real*8 :: varphi_X, X_p, R_port, dist
        real*8 :: varphi_y, ya_plus

        if(ij == 1)then

            ! get initial distribution at age 1 of cash-on-hand
            do iw = 1, NW

                ! get initial cash-on-hand
                X_p = w*eff(1)*zeta(iw)

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                ! get distributional weight
                dist = dist_zeta(iw)

                ! initialize the distribution
                phi_X(1, ixl, 0) = phi_X(1, ixl, 0) + dist*varphi_X
                phi_X(1, ixr, 0) = phi_X(1, ixr, 0) + dist*(1d0-varphi_X)
            enddo

        elseif(ij <= JR-1)then

            ! iterate over yesterdays wealth distribution
            do ia = 0, NA
                do iy = 0, NY

                    do iw = 1, NW
                        do isr = 1, NSR

                            ! get today's cash-on-hand and interpolate
                            R_port = 1d0 + r_f + omega_plus(ij-1, ia, iy)*(mu_r + vtheta(isr))
                            X_p = R_port*a(ia)/eps(isr) + w*eff(ij)*zeta(iw) + ya(iy)
                            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                            ! get future balance of retirement account
                            ya_plus = ya(iy)/eps(isr)

                            ! derive interpolation weights
                            call linint_Grow(ya_plus, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

                            ! get distributional weight
                            dist = dist_zeta(iw)*dist_epsvtheta(isr)

                            ! distribute on today's state space
                            phi_X(ij, ixl, iyl) = phi_X(ij, ixl, iyl) + &
                                        dist*varphi_X*varphi_y*phi_r(ij-1, ia, iy)
                            phi_X(ij, ixr, iyl) = phi_X(ij, ixr, iyl) + &
                                        dist*(1d0-varphi_X)*varphi_y*phi_r(ij-1, ia, iy)
                            phi_X(ij, ixl, iyr) = phi_X(ij, ixl, iyr) + &
                                        dist*varphi_X*(1d0-varphi_y)*phi_r(ij-1, ia, iy)
                            phi_X(ij, ixr, iyr) = phi_X(ij, ixr, iyr) + &
                                        dist*(1d0-varphi_X)*(1d0-varphi_y)*phi_r(ij-1, ia, iy)
                        enddo
                    enddo
                enddo
            enddo

        else

            ! iterate over yesterdays wealth distribution
            do ia = 0, NA
                do iy = 0, NY

                    do isr = 1, NSR

                        ! get today's cash-on-hand
                        R_port = 1d0 + r_f + omega_plus(ij-1, ia, iy)*(mu_r + vtheta(isr))
                        X_p = R_port*a(ia) + pen(ij) + ya(iy)

                        ! derive interpolation weights
                        call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                        ! get distributional weight
                        dist = dist_epsvtheta(isr)

                        ! distribute on today's state space
                        phi_X(ij, ixl, iy) = phi_X(ij, ixl, iy) + dist*varphi_X*phi_r(ij-1, ia, iy)
                        phi_X(ij, ixr, iy) = phi_X(ij, ixr, iy) + dist*(1d0-varphi_X)*phi_r(ij-1, ia, iy)

                    enddo
                enddo
            enddo
        endif

    end subroutine


    ! to calculate distribution after consumption savings decision
    subroutine get_distribution_a(ij)

        implicit none
        integer, intent(in) :: ij
        integer :: ix, iy, ial, iar
        real*8 :: varphi_a

        ! iterate over todays cash on hand and retirement wealth
        do ix = 0, NX
            do iy = 0, NY

                ! interpolate asset decision
                call linint_Grow(a_plus(ij, ix, iy), a_l, a_u, a_grow, NA, ial, iar, varphi_a)

                ! restrict values to grid just in case
                ial = min(ial, NA)
                iar = min(iar, NA)
                varphi_a = min(varphi_a, 1d0)

                ! get distribution over total assets and retirement wealth
                phi_a(ij, ial, iy) = phi_a(ij, ial, iy) + varphi_a*phi_X(ij, ix, iy)
                phi_a(ij, iar, iy) = phi_a(ij, iar, iy) + (1d0-varphi_a)*phi_X(ij, ix, iy)
            enddo
        enddo

    end subroutine


    ! to calculate distribution after annuity savings decision
    subroutine get_distribution_ya(ij)

        implicit none
        integer, intent(in) :: ij
        real*8 :: al_p, ya_p, varphi_a, varphi_y
        integer :: ia, iy, ial, iar, iyl, iyr

        ! iterate over both wealth dimensions
        do ia = 0, NA
            do iy = 0, NY

                ! determine future liquid wealth and pre-return value of ret. account
                al_p = (1d0-omegay_plus(ij, ia, iy))*a(ia)
                ya_p = ya(iy) + omegay_plus(ij, ia, iy)*a(ia)/max(p_a(ij), 1d-10)

                ! derive interpolation weights
                call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
                call linint_Grow(ya_p, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

                ! get distribution over liquid asset and pre-return ret. wealth
                phi_r(ij, ial, iyl) = phi_r(ij, ial, iyl) + varphi_a*varphi_y*phi_a(ij, ia, iy)
                phi_r(ij, iar, iyl) = phi_r(ij, iar, iyl) + (1d0-varphi_a)*varphi_y*phi_a(ij, ia, iy)
                phi_r(ij, ial, iyr) = phi_r(ij, ial, iyr) + varphi_a*(1d0-varphi_y)*phi_a(ij, ia, iy)
                phi_r(ij, iar, iyr) = phi_r(ij, iar, iyr) + (1d0-varphi_a)*(1d0-varphi_y)*phi_a(ij, ia, iy)
            enddo
        enddo

    end subroutine


    ! subroutine for calculating quantities
    subroutine aggregation()

        implicit none
        integer :: ij, ix, ia, iw, iy
        real*8 :: sigma_eta(JJ), mu_exp(JJ), sigma_exp(JJ)
        real*8 :: var_c(JJ), var_y(JJ), var_a(JJ), var_ya(JJ), var_o(JJ)

        ! calculate cohort averages
        y_coh(:)  = 0d0
        c_coh(:)  = 0d0
        ya_coh(:) = 0d0
        a_coh(:)  = 0d0
        o_coh(:)  = 0d0

        do ij = 1, JJ

            do iw = 1, NW
                y_coh(ij) = y_coh(ij) + (w*eff(ij)*zeta(iw) + pen(ij))*dist_zeta(iw)
            enddo

            do iy = 0, NY
                do ix = 0, NX
                    c_coh(ij) = c_coh(ij) + c(ij, ix, iy)*phi_X(ij, ix, iy)
                    ya_coh(ij) = ya_coh(ij)+ ya(iy)*phi_X(ij, ix, iy)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        a_coh(ij)  = a_coh(ij) + a(ia)*phi_r(ij-1, ia, iy)
                        o_coh(ij)  = o_coh(ij) + omega_plus(ij-1, ia, iy)*phi_r(ij-1, ia, iy)
                    enddo
                endif
            enddo
        enddo

        do iy = 0, NY
            do ia = 0, NA
                a_coh(JJ+1)  = a_coh(JJ+1) + a(ia)*phi_r(JJ, ia, iy)
                o_coh(JJ+1)  = o_coh(JJ+1) + omega_plus(JJ, ia, iy)*phi_r(JJ, ia, iy)
            enddo
        enddo

        ! calculate variances
        var_y  = 0d0
        var_c  = 0d0
        var_ya = 0d0
        var_a  = 0d0
        var_o  = 0d0

        do ij = 1, JJ

            do iw = 1, NW
                var_y(ij) = var_y(ij) + (w*eff(ij)*zeta(iw) + pen(ij))**2*dist_zeta(iw)
            enddo

            do iy = 0, NY
                do ix = 0, NX
                    var_c(ij) = var_c(ij) + c(ij, ix, iy)**2*phi_X(ij, ix, iy)
                    var_ya(ij) = var_ya(ij)+ ya(iy)**2*phi_X(ij, ix, iy)
                enddo

                if(ij > 1)then
                    do ia = 0, NA
                        var_a(ij)  = var_a(ij) + a(ia)**2*phi_r(ij-1, ia, iy)
                        var_o(ij)  = var_o(ij) + omega_plus(ij-1, ia, iy)**2*phi_r(ij-1, ia, iy)
                    enddo
                endif
            enddo
        enddo

        var_c  = var_c - c_coh**2
        var_ya = var_ya- ya_coh**2
        var_y  = var_y - y_coh**2
        var_a  = var_a - a_coh(1:JJ)**2
        var_o  = var_o - o_coh(1:JJ)**2

        ! get age dependent variance of eta
        sigma_eta = sigma_eps*(/(dble(min(ij, JR-1)-1), ij=1,JJ)/)

        ! calculate age specific expectations and variance of exp(eta)
        mu_exp = exp(0.5d0*sigma_eta)
        sigma_exp = exp(sigma_eta)*(exp(sigma_eta)-1d0)

        ! add variance effects
        var_y  = mu_exp**2*var_y + sigma_exp*y_coh**2 + sigma_exp*var_y
        var_c  = mu_exp**2*var_c + sigma_exp*c_coh**2 + sigma_exp*var_c
        var_ya = mu_exp**2*var_ya+ sigma_exp*ya_coh**2+ sigma_exp*var_ya
        var_a  = mu_exp**2*var_a + sigma_exp*a_coh(1:JJ)**2 + sigma_exp*var_a


        ! add level effect to averages
        y_coh  = mu_exp*y_coh
        c_coh  = mu_exp*c_coh
        ya_coh = mu_exp*ya_coh
        a_coh(1:JJ)  = mu_exp*a_coh(1:JJ)

        ! calculate coefficients of variation
        cv_y  = sqrt(max(var_y, 0d0))/max(y_coh, 1d-10)
        cv_c  = sqrt(max(var_c, 0d0))/max(c_coh, 1d-10)
        cv_ya = sqrt(max(var_ya,0d0))/max(ya_coh,1d-10)
        cv_a  = sqrt(max(var_a, 0d0))/max(a_coh(1:JJ), 1d-10)
        cv_o  = sqrt(max(var_o, 0d0))/max(o_coh(1:JJ), 1d-10)

    end subroutine


    ! subroutine for writing output
    subroutine output()

        implicit none
        integer :: ij, ixmax(JJ), iamax(JJ), iymax(JJ), ages(JJ+1)

        ! check for the maximium grid points used
        call check_grid_X(ixmax)
        call check_grid_a(iamax)
        call check_grid_ya(iymax)

        ! set up age variable
        ages = (/(20+ij, ij=1,JJ+1)/)

        write(21, '(a,a,a)')' IJ      CONS    INCOME       Y_A    ASSETS', &
            '     OMEGA     CV(C)     CV(Y)    CV(YA)     CV(A)     CV(O)', &
            '     IXMAX     IAMAX     IYMAX'
        do ij = 1, JJ
            write(21,'(i3,10f10.3,3i10)')ages(ij), c_coh(ij), y_coh(ij), ya_coh(ij), a_coh(ij), &
                    o_coh(ij), cv_c(ij), cv_y(ij), cv_ya(ij), cv_a(ij), cv_o(ij), &
                    ixmax(ij), iamax(ij), iymax(ij)
        enddo
        write(21,'(i3,30x,2f10.3)')ages(JJ+1), a_coh(ij), o_coh(ij)
        write(21,'(a/)')'--------------------------------------------------------------------'

        ! plot output
        call plot(dble(ages(1:JJ)), c_coh, legend='Consumption (Mean)')
        call plot(dble(ages(1:JJ)), y_coh, legend='Labor Income (Mean)')
        call plot(dble(ages(1:JJ)), ya_coh, legend='Annuity Payments (Mean)')
        call execplot(xlabel='Age j', ylabel='Consumption/Income/Annuities', ylim=(/0d0, 4d0/))

        call plot(dble(ages(2:JJ+1)), o_coh(2:JJ+1))
        call execplot(xlabel='Age j', ylabel='Portfolio Share', ylim=(/0d0, 1d0/))

        call plot(dble(ages), a_coh, legend='Liquid Assets')
        call plot(dble(ages(1:JJ)), ya_coh, legend='Annuity stream')
        call execplot(xlabel='Age j', ylabel='Assets')

    end subroutine


    ! subroutine that checks for the maximum asset gridpoint used
    subroutine check_grid_a(iamax)

        implicit none
        integer :: iamax(JJ), ij, ia, iy

        iamax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ia = 0, NA
                do iy = 0, NY
                    if(phi_a(ij, ia, iy) > 1d-8)iamax(ij) = ia
                    if(phi_r(ij, ia, iy) > 1d-8)iamax(ij) = ia
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_X(ixmax)

        implicit none
        integer :: ixmax(JJ), ij, ix, iy

        ixmax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do iy = 0, NY
                    if(phi_X(ij, ix, iy) > 1d-8)ixmax(ij) = ix
                enddo
            enddo
        enddo

    end subroutine


    ! subroutine that checks for the maximum gridpoint used
    subroutine check_grid_ya(iymax)

        implicit none
        integer :: iymax(JJ), ij, ix, iy

        iymax = 0
        do ij = 1, JJ

            ! check for the maximum asset grid point used at a certain age
            do ix = 0, NX
                do iy = 0, NY
                    if(phi_X(ij, ix, iy) > 1d-8)iymax(ij) = iy
                enddo
            enddo
        enddo

    end subroutine

end program
