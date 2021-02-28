!##############################################################################
! MODULE globals
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
module globals

    use toolbox

    implicit none

    ! number of years the household retires
    integer, parameter :: JR = 45

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of persistent shock process values
    integer, parameter :: NP = 2

    ! number of transitory shock process values
    integer, parameter :: NS = 7

    ! number of points on the asset grid
    integer, parameter :: NA = 200

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu    = 0.335d0
    real*8, parameter :: beta  = 0.98d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 200d0
    real*8, parameter :: a_grow = 0.05d0

    ! net prices
    real*8 :: r, w

    ! transfer payments (old-age) and survival probabilities
    real*8 :: pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), h_coh(JJ), a_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_l(JJ), cv_h(JJ), corr_hl(JJ)

    ! the permanent shock process
    real*8 :: dist_theta(NP), theta(NP)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 4


    ! demographic and other model parameters
    real*8 :: eff(JJ)
    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, NP, NS)
    real*8 :: c(JJ, 0:NA, NP, NS), l(JJ, 0:NA, NP, NS)
    real*8 :: phi(JJ, 0:NA, NP, NS), V(JJ, 0:NA, NP, NS) = 0d0

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, NP, NS), EV(JJ, 0:NA, NP, NS)
    integer :: ij_com, ia_com, ip_com, is_com
    real*8 :: cons_com, lab_com

contains


    ! the first order condition
    function foc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, a_plus, varphi, tomorrow, wage, available
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus  = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! calculate available resources
        available = (1d0+r)*a(ia_com) + pen(ij_com)

        ! determine labor
        if(ij_com < JR)then
            lab_com = min( max( nu + (1d0-nu)*(a_plus - available)/wage, 0d0), 1d0-1d-10)
        else
            lab_com = 0d0
        endif

        ! calculate consumption
        cons_com = max( (available + wage*lab_com - a_plus) , 1d-10)

        ! calculate linear interpolation for future part of first order condition
        a_plus = max(a_plus, a_l)
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com+1, ial, ip_com, is_com) + &
                  (1d0-varphi)*RHS(ij_com+1, iar, ip_com, is_com)

        ! calculate first order condition for consumption
        foc = margu(cons_com, lab_com)**(-gamma) - tomorrow

    end function


    ! calculates marginal utility of consumption
    function margu(cons, lab)

        implicit none
        real*8, intent(in) :: cons, lab
        real*8 :: margu

        margu = nu*(cons**nu*(1d0-lab)**(1d0-nu))**egam/cons

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, lab, ij, ip, is)

        implicit none
        integer, intent(in) :: ij, ip, is
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
            valuefunc = max(varphi*EV(ij+1, ial, ip, is) + &
                (1d0-varphi)*EV(ij+1, iar, ip, is), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = (c_help**nu*(1d0-l_help)**(1d0-nu))**egam/egam + beta*psi(ij+1)*valuefunc

    end function

end module
