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

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of years the household retires
    integer, parameter :: JR = 45

    ! number of points on the asset grid
    integer, parameter :: NA = 200

    ! number of persistent shock process values
    integer, parameter :: NP = 2

    ! number of transitory shock process values
    integer, parameter :: NS = 7

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: chi = 0.6d0
    real*8, parameter :: echi = 1d0 + 1d0/chi
    real*8, parameter :: nu = 25d0
    real*8, parameter :: beta = 0.98d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! tax progressivity parameters
    real*8, parameter :: tau = 0.2d0
    real*8, parameter :: varrho = 0.8d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 200d0
    real*8, parameter :: a_grow = 0.05d0

    ! the permanent shock process
    real*8 :: dist_theta(NP), theta(NP)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 4

    ! net prices
    real*8 :: r, w

    ! transfer payments (old-age), survival probabilities, productivity
    real*8 :: pen(JJ), psi(JJ+1), psi_cum(JJ+1), eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, NP, NS)
    real*8 :: c(JJ, 0:NA, NP, NS), l(JJ, 0:NA, NP, NS)
    real*8 :: phi(JJ, 0:NA, NP, NS), V_c(JJ, 0:NA, NP, NS) = 0d0, V_l(JJ, 0:NA, NP, NS) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), n_coh(JJ), l_coh(JJ), h_coh(JJ), a_coh(JJ), v_coh(JJ), t_coh(JJ)
    real*8 :: cv_c(JJ), cv_n(JJ), cv_l(JJ), cv_h(JJ), corr_hn(JJ), tax_revenue, labor_supply

    ! endogenous gridpoint variables
    real*8 :: a_endo(0:NA), c_endo(0:NA), l_endo(0:NA)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, NP, NS), EV_c(JJ, 0:NA, NP, NS), EV_l(JJ, 0:NA, NP, NS)
    integer :: ij_com, ia_com, ip_com, is_com
    real*8 :: a_plus_com, lab_com

contains


    ! the first order condition
    function foc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, cons, varphi, tomorrow, wage, expo
        integer :: ial, iar

        ! today's consumption
        cons = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! get exponent
        expo = 1d0/(1d0+1d0/chi-varrho)

        ! get labor supply
        lab_com = max(((1d0-tau)*varrho*wage**varrho/nu)**expo* &
                      max(cons, 1d-10)**(-expo/gamma), 0d0)

        ! calculate tomorrows assets
        a_plus_com = (1d0+r)*a(ia_com) + pen(ij_com) + (1d0-tau)* &
                     (wage*lab_com)**varrho - cons

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(a_plus_com, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com+1, ial, ip_com, is_com) + &
                   (1d0-varphi)*RHS(ij_com+1, iar, ip_com, is_com)

        ! calculate first order condition for consumption
        foc = cons - tomorrow

    end function


    ! the first order condition for liquidity constrained agent
    function foc_liq(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_liq, cons, wage, expo

        ! today's consumption
        cons = x_in

        ! get exponent
        expo = 1d0/(1d0+1d0/chi-varrho)

        ! calculate the wage rate
        wage = w*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! get labor supply
        lab_com = max(((1d0-tau)*varrho*wage**varrho/nu)**expo*max(cons, 1d-10) &
                      **(-expo/gamma), 0d0)

        ! calculate first order condition for consumption
        foc_liq = (1d0+r)*a(ia_com) + pen(ij_com) + (1d0-tau)*(wage*lab_com)**varrho - cons

    end function


    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = max(cons, 1d-10)**(-1d0/gamma)

    end function


    ! calculates the value function of consumption
    function valuefunc_c(a_plus, cons, ij, ip, is)

        implicit none
        integer, intent(in) :: ij, ip, is
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc_c, varphi, c_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc_c = 0d0
        if(ij < JJ)then
            valuefunc_c = max(varphi*EV_c(ij+1, ial, ip, is) + &
                         (1d0-varphi)*EV_c(ij+1, iar, ip, is), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc_c = c_help**egam/egam + beta*psi(ij+1)*valuefunc_c

    end function


    ! calculates the value function of labor
    function valuefunc_l(a_plus, lab, ij, ip, is)

        implicit none
        integer, intent(in) :: ij, ip, is
        real*8, intent(in) :: a_plus, lab
        real*8 :: valuefunc_l, varphi, l_help
        integer :: ial, iar

        ! check whether labor is greater than 0
        l_help = max(lab, 0d0)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc_l = 0d0
        if(ij < JJ)then
            valuefunc_l = max(varphi*EV_l(ij+1, ial, ip, is) + &
                         (1d0-varphi)*EV_l(ij+1, iar, ip, is), 1d-10)**echi/echi
        endif

        ! add todays part and discount
        valuefunc_l = nu*l_help**echi/echi + beta*psi(ij+1)*valuefunc_l

    end function

end module
