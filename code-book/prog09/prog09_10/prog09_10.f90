!##############################################################################
! PROGRAM HetAg_Gov
!
! ## Heterogeneous agents and government activity (optimal quantity of debt)
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
include "prog09_10m.f90"

program HetAg_Gov

    use globals
    use toolbox

    implicit none
    integer, parameter :: NB = 4
    integer :: ib
    real*8 :: blevels(0:NB), Wlevels(0:NB), W0

    ! initialize variables
    call initialize()

    ! benchmark equilibrium
    by = 0.6d0

    ! set initial guess of the interest rate
    r = 0.042d0

    ! solve initial equilibrium
    call solve_model()

    ! store initial uility level
    W0 = sum(phi(:,:)*V(:,:))

    ! set up number of debt levels
    call grid_Cons_Equi(blevels, 0.5d0, 1.5d0)

    ! calculate equilibrium with different debt levels
    do ib = 0, NB
        by = blevels(ib)
        call solve_model()
        Wlevels(ib) = sum(phi(:,:)*V(:,:))
    enddo

    ! plot equivalent variation
    call plot(blevels, ((Wlevels/W0)**(1d0/egam)-1d0)*100d0)
    call execplot

contains


    ! For solving the model
    subroutine solve_model

        use toolbox

        implicit none
        real*8 :: x_in
        logical :: check

        ! start timer
        call tic()

        write(*, '(a)')'       K/Y         r     tau_r      diff_K'

        ! set tolerance level for outer optimization
        call settol_root(sig_out)

        ! set initial guess
        x_in = r
        check = .false.

        ! solve asset market equilibrium
        call fzero(x_in, asset_market, check)

        if(check)write(*,'(a)')'ERROR IN ROOTFINDING'

        call toc

        ! write output to the screen
        call output

    end subroutine


    ! For initializing variables
    subroutine initialize()

        use toolbox

        implicit none
        integer :: is

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi, weights)
        eta = exp(eta)

        ! calculate aggregate labor supply
        LL = sum(eta*weights)

        ! initialize grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize policy function
        do is = 1, ns
            c(:, is) = 0.04d0*a(:) + eta(is)
            c_new(:, is) = c(:, is)
        enddo

        ! get an initial guess for the distribution function
        phi = 1d0/dble((NA+1)*NS)

    end subroutine


    ! For calculating the value function
    subroutine value_function()

        use toolbox

        implicit none
        integer :: iter, ia, is, is_p
        real*8 :: V_new(0:NA, NS), aplus, con_lev
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
                aplus = max(((1d0+rn)*a(ia)+wn*eta(is)-c(ia,is))/(1d0+growth), a_l)

                ! determine the gridpoints in between this decision lies and share
                call linint_Grow(aplus, a_l, a_u, a_grow, NA, &
                    ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo

        ! initialize value function
        V = c

        ! iterate until the value function converges
        do iter = 1, itermax

            ! calculate new value function
            do ia = 0, NA
                do is = 1, NS

                    ! interpolate over all future states
                    V_new(ia, is) = 0d0
                    do is_p = 1, NS
                        V_new(ia, is) = V_new(ia, is) + pi(is, is_p)* &
                            (varphi(ia,is)*V(ial(ia,is),is_p) + &
                            (1d0-varphi(ia,is))*V(iar(ia,is),is_p))**egam
                    enddo

                    ! add current utility
                    V_new(ia, is) = (c(ia,is)**egam + beta*V_new(ia, is))**(1d0/egam)
                enddo
            enddo

            con_lev = maxval(abs(V_new(:, :) - V(:, :))/max(abs(V(:, :)), 1d-10))

            ! update value function
            V = V_new

            ! check for convergence
            if(con_lev < sig_in)then
                V = V**egam/egam
                return
            endif
        enddo

    end subroutine


    ! For creating output plots
    subroutine output()

        implicit none
        integer :: ia, is
        real*8 :: Std_c, Std_a, Std_eta, util

        ! aggregate consumption and investment
        CC = sum(phi*c)
        II = (growth+delta)*KK

        call value_function()
        util = sum(phi(:,:)*V(:,:))

        ! calculate coefficients of variation
        Std_eta = 0d0
        Std_c   = 0d0
        Std_a   = 0d0
        do ia = 0, NA
            do is = 1, NS
                Std_eta = Std_eta + phi(ia, is)*(eta(is)   - LL)**2
                Std_c   = Std_c   + phi(ia, is)*(c(ia, is) - CC)**2
                Std_a   = Std_a   + phi(ia, is)*(a(ia)     - KK)**2
            enddo
        enddo
        Std_eta = sqrt(Std_eta)/LL
        Std_c   = sqrt(Std_c)  /CC
        Std_a   = sqrt(Std_a)  /AA

        ! write output to screen
        write(*, '(//a)')'AIYAGARI ECONOMY OUTPUT'
        write(*,'(/a)')'PROD :       K/Y       L/Y         r         w'
        write(*,'(6x,4f10.4)')(/KK, LL/)/YY*100d0, r*100d0, w
        write(*,'(/a)')'GOODS:       C/Y       I/Y       G/Y       D/Y'
        write(*,'(6x,4f10.4)')(/CC, II, GG, CC+II+GG/)/YY*100d0
        write(*,'(/a)')'GOVMT:       B/Y       G/Y     tau_w     tau_r'
        write(*,'(6x,5f10.4)')(/BB, GG/)/YY*100d0, tau_w, tau_r
        write(*,'(/a)')'               c         a       eta     util'
        write(*,'(a,4f10.4)')'LEV  :', CC, AA, LL, util
        write(*,'(a,3f10.4)')'STD  :', Std_c, Std_a, Std_eta
        write(*,'(//)')

    end subroutine

end program
