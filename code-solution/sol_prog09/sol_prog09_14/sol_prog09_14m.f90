!##############################################################################
! MODULE globals
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
module globals

    use toolbox

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: nu = 0.3d0
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 0.08d0
    real*8, parameter :: rho = 0.6d0
    real*8, parameter :: sigma_eps = 0.04d0*(1d0-rho**2)

    ! numerical parameters
    real*8, parameter :: a_l = 0d0
    real*8, parameter :: a_u = 100d0
    real*8, parameter :: a_grow = 0.01d0
    real*8, parameter :: sig_in = 1d-8
    real*8, parameter :: sig_out = 1d-6
    integer, parameter :: itermax = 50000

    ! macroeconomic variables
    real*8 :: r, w, KK, AA, LL, YY, CC, II

    ! the shock process
    integer, parameter :: NS = 7
    real*8 :: pi(NS, NS), eta(NS), weights(NS)

    ! policy function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA, NS), l(0:NA, NS), leff(0:NA, NS), aplus(0:NA, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NA, NS), coeff_c(NA+3, NS), coeff_l(NA+3, NS)
    real*8 :: RHS(0:NA, NS), coeff_r(NA+3, NS)
    logical :: check

    ! variables to communicate with function
    real*8 :: a_com, l_com, aplus_com
    integer :: is_com

    ! variables to numerically determine the steady state distribution
    real*8 :: phi(0:NA, NS), phi_new(0:NA, NS)

contains


    ! the function to determine equilibrium capital in the aiyagari model
    function asset_market(x)

        real*8, intent(in) :: x
        real*8 :: asset_market, KL

        ! set capital to labor ratio
        KL = x

        ! calculate factor prices
        w = (1d0-alpha)*(KL)**alpha
        r = alpha*(1d0/KL)**(1d0-alpha) - delta

        ! solve the household value function iteration problem
        call solve_household()

        ! calculate the invariant distribution directly
        call get_distribution()

        ! calculate macroeconomic variables
        AA = sum(a*sum(phi, 2))
        LL = sum(leff*phi)
        KK = KL*LL
        YY = KK**alpha*LL**(1d0-alpha)

        ! write output to screen
        write(*,'(2f10.4,f12.8)')KK/YY*100d0, r*100d0, AA - KK

        ! check for difference between input and computed capital
        asset_market = (AA - KK)/AA

    end function


    ! determines the solution to the household optimization problem
    subroutine solve_household()

        implicit none
        integer :: ia, is, iter
        real*8 :: con_lev, p

        ! set tolerance level for interior optimization
        call settol_root(1d-12)

        ! do a policy function iteration with rootfinding
        do iter = 1, itermax

            ! interpolate coefficients
            call interpolate()

            ! solve the household problem for all gridpoints
            do ia = 0, NA
                do is = 1, NS

                    ! initialize starting value and communicate resource level
                    p = c(ia, is)
                    a_com = a(ia)
                    is_com = is

                    ! find the optimal consumption level
                    call fzero(p, foc, check)
                    if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

                    ! get optimal consumption function
                    c_new(ia, is) = p
                    l(ia, is) = l_com
                    aplus(ia, is) = aplus_com

                    ! check for borrowing constraint
                    if(aplus(ia, is) < 0d0)then
                        aplus(ia, is) = 0d0
                        l(ia, is) = max(nu-(1d0-nu)*(1d0+r)*a_com/(w*eta(is_com)), 0d0)
                        c_new(ia, is) = max((1d0+r)*a_com + w*eta(is_com)*l(ia, is), 1d-10)
                    endif

                    leff(ia, is) = eta(is)*l(ia, is)

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
        real*8 :: foc, c_act, aplus, varphi
        integer :: ial, iar

        ! calculate consumption level
        c_act = max(x_in, 1d-10)

        ! calculate optimal labor supply
        l_com = max(1d0-(1d0-nu)/nu*c_act/(w*eta(is_com)), 0d0)

        ! future assets
        aplus_com = (1d0+r)*a_com + w*eta(is_com)*l_com - c_act

        ! calculate future expected marginal utility
        call linint_Grow(aplus_com, a_l, a_u, a_grow, NA, ial, iar, varphi)
        foc = varphi*RHS(ial, is_com) + (1d0-varphi)*RHS(iar, is_com)

        ! evaluate first order condition
        foc = (nu*(c_act**nu*(1d0-l_com)**(1d0-nu))**egam/c_act)**(-gamma) - foc

    end function


    ! For interpolating the RHS of the first order condition
    subroutine interpolate()

        implicit none

        integer :: is, ia, is_p
        real*8 :: marg_util

        do is = 1, NS
            do ia = 0, NA

                ! calculate the RHS of the first order condition
                RHS(ia, is) = 0d0
                do is_p = 1, NS
                    marg_util = nu*(c(ia, is_p)**nu*(1d0-l(ia, is_p))**(1d0-nu))**egam/c(ia, is_p)
                    RHS(ia, is) = RHS(ia, is) + pi(is, is_p)*marg_util
                enddo
                RHS(ia, is) = (beta*(1d0+r)*RHS(ia, is))**(-gamma)
            enddo

            ! interpolate
            call spline_interp(RHS(:, is), coeff_r(:, is))
        enddo

    end subroutine


    ! calculates the invariant distribution of households
    subroutine get_distribution()

        use toolbox

        implicit none
        integer :: ia, is, iter, ial, iar, is_p
        real*8 :: varphi, con_lev

        ! iterate until the distribution function converges
        do iter = 1, itermax

            phi_new = 0d0

            do ia = 0, NA
                do is = 1, NS

                    ! calculate where this guy would go
                    call linint_Grow(aplus(ia, is), a_l, a_u, a_grow, NA, ial, iar, varphi)
                    varphi = max(min(varphi, 1d0), 0d0)

                    do is_p = 1, NS
                        phi_new(ial, is_p) = phi_new(ial, is_p) + &
                            pi(is, is_p)*varphi*phi(ia, is)
                        phi_new(iar, is_p) = phi_new(iar, is_p) + &
                            pi(is, is_p)*(1d0-varphi)*phi(ia, is)
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

        write(*, *)'Distribution Function Iteration did not converge'

    end subroutine

end module
