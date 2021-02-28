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

    ! number of years the household retires
    integer, parameter :: JR = 45

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 7

    ! number of transitory (epsilon) shocks
    integer, parameter :: NS = 7

    ! number of rate of return (vtheta) shocks
    integer, parameter :: NR = 7

    ! number of eps-vtheta shocks
    integer, parameter :: NSR = NS*NR

    ! number of eta shocks
    integer, parameter :: NE = 1000

    ! number of points on the asset grid
    integer, parameter :: NA = 200

    ! number of points on the cash-on-hand-grid
    integer, parameter :: NX = 200

    ! household preference parameters
    real*8, parameter :: gamma = 0.10d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu = 0.4d0
    real*8, parameter :: beta = 0.96d0

    ! risk processes
    real*8, parameter :: sigma_zeta   = 0.0738d0
    real*8, parameter :: sigma_eps    = 0.0106d0
    real*8, parameter :: sigma_vtheta = 0.157d0**2d0
    real*8, parameter :: rho          = 0.00d0

    ! risk free rate and risk premium
    real*8, parameter :: r_f  = 0.02d0
    real*8, parameter :: mu_r = 0.04d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 200d0
    real*8, parameter :: a_grow = 0.04d0

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.04d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! should cohort averages and variance be calculated analytically
    logical, parameter :: analytical = .true.

    ! lower and upper bound for X-grid
    real*8 :: X_l, X_u

    ! lower and upper bound for the eta-grid
    real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w, eff(JJ), pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), l_coh(JJ), y_coh(JJ), a_coh(JJ), o_coh(JJ)
    real*8 :: cv_c(JJ), cv_l(JJ), cv_y(JJ), cv_a(JJ), cv_o(JJ)

    ! individual variables
    real*8 :: X(0:NX), a(0:NA)
    real*8 :: c(JJ, 0:NX, NW), l(JJ, 0:NX, NW), a_plus(JJ, 0:NX, NW), V(JJ, 0:NX, NW) = 0d0
    real*8 :: omega_plus(JJ, 0:NA), Q(JJ, 0:NA)
    real*8 :: phi_X(JJ, 0:NX, NW), phi_a(JJ, 0:NA)
    real*8 :: eta(JJ, 0:NE), phi_e(JJ, 0:NE)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA)
    integer :: ij_com, ix_com, ia_com, iw_com
    real*8 :: cons_com, lab_com

contains


    ! the first order condition regarding consumption
    function foc_cons(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_cons, a_plus, varphi, wage, tomorrow
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*zeta(iw_com)

        ! get labor supply
        if(ij_com < JR)then
            lab_com = min(max(nu + (1d0-nu)/wage*(a_plus - X(ix_com)), 0d0), 1d0-1d-10)
        else
            lab_com = 0d0
        endif

        ! calculate consumption
        cons_com = max(X(ix_com) + wage*lab_com - a_plus, 1d-10)

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com, ial) + (1d0-varphi)*RHS(ij_com, iar)

        ! calculate first order condition for consumption
        foc_cons = margu(cons_com, lab_com)**(-gamma) - tomorrow

    end function


    ! the first order condition regarding portfolio choice
    function foc_port(p)

        implicit none
        real*8, intent(in) :: p
        real*8 :: foc_port, omega_p, R_port, X_p, c_p, l_p, varphi, dist
        integer :: ixl, ixr, iw, isr

        ! store portfolio share
        omega_p = p

        foc_port = 0d0
        if(ij_com+1 >= JR)then
            do isr = 1, NSR

                ! get return on the portfolio
                R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))

                ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                X_p = R_port*a(ia_com) + pen(ij_com+1)

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! calculate consumption and FOC
                c_p = varphi      *c(ij_com+1, ixl, 1) + &
                      (1d0-varphi)*c(ij_com+1, ixr, 1)
                c_p = max(c_p, 1d-10)
                foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*margu(c_p, 1d-10)
            enddo
        else
            do iw = 1, NW
                do isr = 1, NSR

                    ! get return on the portfolio
                    R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))

                    ! get tomorrow's cash on hand
                    X_p = R_port*a(ia_com)/eps(isr)

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_epsvtheta(isr)

                    ! calculate consumption
                    c_p = varphi      *c(ij_com+1, ixl, iw) + &
                          (1d0-varphi)*c(ij_com+1, ixr, iw)
                    c_p = max(c_p, 1d-10)

                    ! get labor supply and FOC
                    l_p = varphi      *l(ij_com+1, ixl, iw) + &
                          (1d0-varphi)*l(ij_com+1, ixr, iw)
                    l_p = min(max(l_p, 0d0), 1d0-1d-10)
                    foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*margu(eps(isr)*c_p, l_p)
                enddo
            enddo
        endif

    end function


    ! calculates marginal utility of consumption
    function margu(cons, lab)

        implicit none
        real*8, intent(in) :: cons, lab
        real*8 :: margu

        margu = nu*(cons**nu*(1d0-lab)**(1d0-nu))**egam/cons

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, lab, ij)

        implicit none
        integer, intent(in) :: ij
        real*8, intent(in) :: a_plus, cons, lab
        real*8 :: valuefunc, varphi, c_help, l_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)
        l_help = min(max(lab, 0d0), 1d0-1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*Q(ij, ial) + (1d0-varphi)*Q(ij, iar), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = (c_help**nu*(1d0-l_help)**(1d0-nu))**egam/egam + beta*psi(ij+1)*valuefunc

    end function

end module
