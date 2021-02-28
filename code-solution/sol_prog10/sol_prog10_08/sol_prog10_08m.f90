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

    ! year the household retires
    integer, parameter :: JR = 45

    ! number of points on the asset grid
    integer, parameter :: NA = 200

    ! number of health states
    integer, parameter :: NM = 1

    ! number of transitory shock process values
    integer, parameter :: NS = 7

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu = 0.335d0
    real*8, parameter :: beta = 0.98d0

    ! health cost parameters
    real*8, parameter :: phi_l = 0d0
    real*8, parameter :: chi = 0.8d0
    real*8, parameter :: varrho = 0.2d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 200d0
    real*8, parameter :: a_grow = 0.05d0

    ! probability of bad health
    real*8 :: dist_m(0:NM), pi_m(JJ, 0:NM, 0:NM)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 4

    ! net prices
    real*8 :: r, w

    ! transfer payments (old-age), survival probabilities, productivity
    real*8 :: pen(JJ), psi(JJ+1, 0:NM), eff(JJ), varrho_m(0:NM)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, 0:NM, NS), c(JJ, 0:NA, 0:NM, NS), l(JJ, 0:NA, 0:NM, NS)
    real*8 :: phi(JJ, 0:NA, 0:NM, NS), V(JJ, 0:NA, 0:NM, NS) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ, 0:NM+1), l_coh(JJ, 0:NM+1), y_coh(JJ, 0:NM+1), a_coh(JJ, 0:NM+1), v_coh(JJ, 0:NM+1)
    real*8 :: cv_c(JJ, 0:NM+1), cv_y(JJ, 0:NM+1), frac_phi(JJ, 0:NM+1)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, 0:NM, NS), EV(JJ, 0:NA, 0:NM, NS)
    integer :: ij_com, ia_com, im_com, is_com
    real*8 :: cons_com, lab_com

contains


    ! the first order condition
    function foc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, a_plus, varphi, tomorrow, wage, available
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*eta(is_com)*varrho_m(im_com)

        ! calculate available resources
        available = (1d0+r)*a(ia_com) + pen(ij_com)

        ! determine labor
        if(ij_com < JR)then
            lab_com = min(max((1d0-nu)*(a_plus-available)/wage + &
                               nu*(1d0-phi_l*dble(im_com)), 1d-10) &
                               , 1d0-phi_l*dble(im_com)-1d-10)
        else
            lab_com = 0d0
        endif

        ! calculate consumption
        cons_com = max(available + wage*lab_com - a_plus, 1d-10)

        ! calculate linear interpolation for future part of first order condition
        a_plus = max(a_plus, a_l)
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com+1, ial, im_com, is_com) + &
                            (1d0-varphi)*RHS(ij_com+1, iar, im_com, is_com)

        ! calculate first order condition for consumption
        foc = margu(cons_com, lab_com, im_com)**(-gamma) - tomorrow

    end function


    ! calculates marginal utility of consumption
    function margu(cons, lab, im)

        implicit none
        real*8, intent(in) :: cons, lab
        integer, intent(in) :: im
        real*8 :: margu

        margu = nu*(cons**nu*(1d0-lab-phi_l*dble(im))**(1d0-nu))**egam/cons

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, lab, ij, im, is)

        implicit none
        integer, intent(in) :: ij, im, is
        real*8, intent(in) :: a_plus, cons, lab
        real*8 :: valuefunc, varphi, c_help, l_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)
        l_help = min(max(lab, 0d0), 1d0-phi_l*dble(im)-1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*EV(ij+1, ial, im, is) + &
                            (1d0-varphi)*EV(ij+1, iar, im, is), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = (c_help**nu*(1d0-l_help-phi_l*dble(im))**(1d0-nu))**egam/egam + beta*psi(ij+1, im)*valuefunc

    end function

end module
