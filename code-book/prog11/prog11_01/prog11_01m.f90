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
    integer, parameter :: JJ = 12

    ! number of years the household retires
    integer, parameter :: JR = 10

    ! number of persistent shock process values
    integer, parameter :: NP = 2

    ! number of transitory shock process values
    integer, parameter :: NS = 5

    ! number of points on the asset grid
    integer, parameter :: NA = 100

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu    = 0.335d0
    real*8, parameter :: beta  = 0.998**5

    ! household risk process
    real*8, parameter :: sigma_theta = 0.23d0
    real*8, parameter :: sigma_eps   = 0.05d0
    real*8, parameter :: rho         = 0.98d0

    ! production parameters
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 1d0-(1d0-0.0823d0)**5
    real*8, parameter :: Omega = 1.60d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 35d0
    real*8, parameter :: a_grow = 0.05d0

    ! demographic parameters
    real*8, parameter :: n_p   = (1d0+0.01d0)**5-1d0

    ! simulation parameters
    real*8, parameter :: damp    = 0.30d0
    real*8, parameter :: sig     = 1d-4
    integer, parameter :: itermax = 50

    ! counter variables
    integer :: iter

    ! macroeconomic variables
    real*8 :: r, rn, w, wn, p
    real*8 :: KK, AA, BB, LL, HH
    real*8 :: YY, CC, II, GG, INC

    ! government variables
    real*8 :: tauc, tauw, taur, taup, kappa
    real*8 :: gy, by, pen(JJ), PP, taxrev(4)
    integer :: tax
    logical :: reform_on

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), a_coh(JJ), v_coh(JJ)

    ! the shock process
    real*8 :: dist_theta(NP), theta(NP)
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 3

    ! demographic and other model parameters
    real*8 :: m(JJ), eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, NP, NS)
    real*8 :: c(JJ, 0:NA, NP, NS), l(JJ, 0:NA, NP, NS)
    real*8 :: phi(JJ, 0:NA, NP, NS), V(JJ, 0:NA, NP, NS) = 0d0

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, NP, NS), ev(JJ, 0:NA, NP, NS)
    integer :: ij_com, ia_com, ip_com, is_com, it_com
    real*8 :: cons_com, lab_com, DIFF, INC_init

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
        wage = wn*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! calculate available resources
        available = (1d0+rn)*a(ia_com) + pen(ij_com)

        ! determine labor
        if(ij_com < JR)then
            lab_com = min( max( nu + (1d0-nu)*(a_plus - available)/wage, 0d0), 1d0-1d-10)
        else
            lab_com = 0d0
        endif

        ! calculate consumption
        cons_com = max( (available + wage*lab_com - a_plus)/p , 1d-10)

        ! calculate linear interpolation for future part of first order condition
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

        margu = nu/p*(cons**nu*(1d0-lab)**(1d0-nu))**egam/cons

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
        l_help = min(max(lab, 0d0),1d0-1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*EV(ij+1, ial, ip, is) + &
                       (1d0-varphi)*EV(ij+1, iar, ip, is), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = (c_help**nu*(1d0-l_help)**(1d0-nu))**egam/egam + beta*valuefunc

    end function

end module
