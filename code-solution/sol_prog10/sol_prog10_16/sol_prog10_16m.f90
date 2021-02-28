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

    ! number of white noise shocks
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
    real*8, parameter :: a_u    = 500d0
    real*8, parameter :: a_grow = 0.04d0

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.04d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! should cohort averages and variance be calculated analytically
    logical, parameter :: analytical = .false.

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
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ), o_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_a(JJ), cv_o(JJ)

    ! individual variables
    real*8 :: X(0:NX), a(0:NA)
    real*8 :: c(JJ, 0:NX), a_plus(JJ, 0:NX), V(JJ, 0:NX) = 0d0
    real*8 :: omega_plus(JJ, 0:NA), Q(JJ, 0:NA)
    real*8 :: phi_X(JJ, 0:NX), phi_a(JJ, 0:NA)
    real*8 :: eta(JJ, 0:NE), phi_e(JJ, 0:NE)

    ! numerical variables
    integer :: ij_com, ix_com, ia_com
    real*8 :: cons_com

contains


    ! the value fucntion regarding consumption
    function value_cons(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: value_cons, a_plus, cons_lim, varphi, tomorrow
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate consumption
        cons_com = max(X(ix_com) - a_plus, 1d-10)
        cons_lim = max(cons_com, 1d-10)

        ! calculate linear interpolation for future value function
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*Q(ij_com, ial) + (1d0-varphi)*Q(ij_com, iar)

        ! calculate value function for consumption
        value_cons = -(cons_com**egam/egam + beta*psi(ij_com+1)*tomorrow**egam/egam)
        value_cons = value_cons + 100d0*abs(cons_com-cons_lim)

    end function


    ! the value function regarding portfolio choice
    function value_port(p)

        implicit none
        real*8, intent(in) :: p
        real*8 :: value_port, omega_p, R_port, X_p, income, V_p, varphi, dist
        integer :: ixl, ixr, iw, isr

        ! store portfolio share
        omega_p = p

        value_port = 0d0
        if(ij_com+1 >= JR)then
            do isr = 1, NSR

                ! get return on the portfolio
                R_port = 1d0 + r_f + omega_p*(mu_r+vtheta(isr))

                ! get tomorrow's cash-on-hand (epsilon^+ = 0)
                X_p = R_port*a(ia_com) + pen(ij_com+1)

                ! derive interpolation weights
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! calculate the future value function
                V_p = varphi      *V(ij_com+1, ixl) + &
                      (1d0-varphi)*V(ij_com+1, ixr)
                V_p = max(V_p, 1d-10)
                value_port = value_port + dist*V_p**egam/egam
            enddo
        else
            do iw = 1, NW
                do isr = 1, NSR

                    ! get return on the portfolio
                    R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))

                    ! derive labor earnings
                    income = w*eff(ij_com+1)*zeta(iw)

                    ! get tomorrow's cash on hand
                    X_p = R_port*a(ia_com)/eps(isr) + income

                    ! derive interpolation weights
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_epsvtheta(isr)

                    ! calculate the future value function
                    V_p = varphi      *V(ij_com+1, ixl) + &
                          (1d0-varphi)*V(ij_com+1, ixr)
                    V_p = max(V_p, 1d-10)
                    value_port = value_port + dist*(eps(isr)*V_p)**egam/egam

                enddo
            enddo
        endif
        value_port = -value_port

    end function

end module
