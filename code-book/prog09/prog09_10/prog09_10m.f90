!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 1d0/1.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: growth = 0.0185d0
    real*8, parameter :: beta = 0.98d0*(1d0+growth)**egam
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: delta = 0.075d0
    real*8, parameter :: rho = 0.6d0
    real*8, parameter :: sigma_eps = 0.09d0*(1d0-rho**2)
    real*8, parameter :: tau_w = 0.27d0
    real*8, parameter :: gy    = 0.217d0

    ! numerical parameters
    real*8, parameter :: a_l = 0d0
    real*8, parameter :: a_u = 100d0
    real*8, parameter :: a_grow = 0.01d0
    real*8, parameter :: sig_in = 1d-10
    real*8, parameter :: sig_out = 1d-6
    integer, parameter :: itermax = 50000

    ! macroeconomic variables
    real*8 :: KK, AA, BB, LL, YY, CC, II, GG
    real*8 :: r, w, rn, wn, tau_r, by

    ! the shock process
    integer, parameter :: NS = 7
    real*8 :: pi(NS, NS), eta(NS), weights(NS)

    ! policy and value function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA, NS), V(0:NA, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NA, NS), RHS(0:NA, NS)
    logical :: check

    ! variables to communicate with function
    real*8 :: a_com
    integer :: is_com

    ! variables to numerically determine the steady state distribution
    real*8 :: phi(0:NA, NS), phi_new(0:NA, NS)

contains


    ! calculate the difference on the asset market given an interest rate
    function asset_market(r_input)

        implicit none
        real*8, intent(in) :: r_input
        real*8 :: asset_market

        ! set the interest rate
        r = r_input

        ! calculate the corresponding capital stock
        KK = (alpha/(r+delta))**(1d0/(1d0-alpha))*LL

        ! get wages and output
        w = (1d0-alpha)*(KK/LL)**alpha
        YY = KK**alpha*LL**(1d0-alpha)

        ! calculate government quantities
        GG = gy*YY
        BB = by*YY

        ! back out aggregate asset holdings
        AA = KK + BB

        ! get budget balancing tax rates
        tau_r = (GG + (r-growth)*BB - tau_w*w*LL)/(r*AA)

        ! define net prices
        rn = r*(1d0-tau_r)
        wn = w*(1d0-tau_w)

        ! solve the household value function iteration problem
        call solve_household()

        ! calculate the invariant distribution directly
        call get_distribution()

        ! get aggregate assets
        AA = sum(a*sum(phi, 2))

        ! get the difference between asset demand and supply
        asset_market = (AA - KK - BB)/AA

        write(*,'(3f10.4,f12.8)')KK/YY*100d0, r*100d0, tau_r*100d0, AA-KK-BB

    end function


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        use toolbox

        implicit none
        integer :: ia, is, iter
        real*8 :: con_lev, x_in

        ! set tolerance level for interior optimization
        call settol_root(sig_in)

        ! do a policy function iteration with rootfinding
        do iter = 1, itermax

            ! interpolate coefficients
            call interpolate()

            ! solve the household problem for all gridpoints
            do ia = 0, NA
                do is = 1, NS

                    ! initialize starting value/communication variables
                    x_in = c(ia, is)
                    a_com = a(ia)
                    is_com = is

                    ! find the optimal consumption level
                    call fzero(x_in, foc, check)
                    if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                    ! check for borrowing constraint
                    if(x_in > (1d0+rn)*a_com+wn*eta(is_com) - a_l) &
                       x_in = (1d0+rn)*a_com+wn*eta(is_com) - a_l

                    ! get optimal consumption and value function
                    c_new(ia, is) = x_in

                enddo
            enddo

            ! get convergence level
            con_lev = maxval(abs(c_new(:, :) - c(:, :))/max(abs(c(:, :)), 1d-10))

            ! check for convergence
            if(con_lev < sig_in)then
                return
            endif

            c = c_new
        enddo

        write(*,*)'Policy Function Iteration did not converge'

    end subroutine


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, aplus, varphi
        integer :: ial, iar

        ! future assets
        aplus = ((1d0+rn)*a_com + wn*eta(is_com) - x_in)/(1d0+growth)

        ! calculate future expected marginal utility
        call linint_Grow(aplus, a_l, a_u, a_grow, NA, ial, iar, varphi)
        foc = varphi*RHS(ial,is_com) + (1d0-varphi)*RHS(iar,is_com)

        ! get first order condition
        foc = x_in - foc

    end function


    ! Set up data for interpolation
    subroutine interpolate()

        implicit none
        integer :: is, ia, is_p

        do is = 1, NS
            do ia = 0, NA

                ! calculate the RHS of the first order condition
                RHS(ia, is) = 0d0
                do is_p = 1, NS
                    RHS(ia, is) = RHS(ia, is) + pi(is, is_p)*c_new(ia, is_p)**(-1d0/gamma)
                enddo
                RHS(ia, is) = (beta*(1d0+rn)*RHS(ia, is))**(-gamma)
            enddo
        enddo

    end subroutine


    ! calculates the invariant distribution of households
    subroutine get_distribution()

        use toolbox

        implicit none
        integer :: ia, is, iter, is_p
        real*8 :: aplus, con_lev
        integer :: ial(0:NA, NS), iar(0:NA, NS)
        real*8 :: varphi(0:NA, NS)

        ! get the interpolation shares and points
        do ia = 0, NA
            do is = 1, NS

                ! calculate where this guy would go
                aplus = ((1d0+rn)*a(ia) + wn*eta(is) - c(ia,is))/(1d0+growth)
                aplus = min(max(aplus, a_l), a_u)

                ! determine the gridpoints in between this decision lies
                call linint_Grow(aplus, a_l, a_u, a_grow, NA, &
                    ial(ia, is), iar(ia, is), varphi(ia, is))
            enddo
        enddo

        ! iterate until the distribution function converges
        do iter = 1, itermax

            phi_new = 0d0

            do ia = 0, NA
                do is = 1, NS

                    do is_p = 1, NS
                        phi_new(ial(ia,is), is_p) = phi_new(ial(ia,is), is_p) + &
                            pi(is,is_p)*varphi(ia,is)*phi(ia,is)
                        phi_new(iar(ia,is), is_p) = phi_new(iar(ia,is), is_p) + &
                            pi(is,is_p)*(1d0-varphi(ia,is))*phi(ia,is)
                    enddo
                enddo
            enddo

            con_lev = maxval(abs(phi_new(:, :) - phi(:, :))/max(abs(phi(:, :)), 1d-10))

            ! update distribution
            phi = phi_new

            ! check for convergence
            if(con_lev < sig_in)then
                return
            endif
        enddo

        write(*,*)'Distribution Function Iteration did not converge'

    end subroutine

end module
