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

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! year the household retires
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
    real*8, parameter :: beta  = 0.98d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 600d0
    real*8, parameter :: a_grow = 0.05d0

    ! the permanent shock process
    real*8 :: dist_theta(NP), theta(NP)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 4

    ! net prices
    real*8 :: r, w

    ! transfer payments (old-age), survival probabilities, productivity
    real*8 :: pen(JJ), psi(JJ+1), eff(JJ), a_bor(JJ, NP)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, NP, NS), c(JJ, 0:NA, NP, NS)
    real*8 :: phi(JJ, 0:NA, NP, NS), V(JJ, 0:NA, NP, NS) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, NP, NS), EV(JJ, 0:NA, NP, NS)
    integer :: ij_com, ia_com, ip_com, is_com
    real*8 :: cons_com

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
        available = (1d0+r)*(a(ia_com) + a_bor(ij_com, ip_com)) + wage + pen(ij_com)

        ! calculate consumption
        cons_com = available - (a_plus + a_bor(ij_com+1, ip_com))

        ! calculate linear interpolation for future part of foc
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com+1, ial, ip_com, is_com) + &
                            (1d0-varphi)*RHS(ij_com+1, iar, ip_com, is_com)

        ! calculate first order condition for consumption
        foc = cons_com - tomorrow

    end function


    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = cons**(-1d0/gamma)

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, ij, ip, is)

        implicit none
        integer, intent(in) :: ij, ip, is
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc, varphi, c_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*EV(ij+1, ial, ip, is) + &
                       (1d0-varphi)*EV(ij+1, iar, ip, is), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = c_help**egam/egam + beta*psi(ij+1)*valuefunc

    end function

end module
