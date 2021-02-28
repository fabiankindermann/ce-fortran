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
    integer, parameter :: NS = 5

    ! number of points on the asset grid
    integer, parameter :: NA = 80

    ! number of points on the human capital grid
    integer, parameter :: NH = 60

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam  = 1d0 - 1d0/gamma
    real*8, parameter :: nu    = 0.12d0
    real*8, parameter :: beta  = 0.98d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

    ! parameters governing the human capital process of women
    real*8, parameter :: xi(2) = (/0.05312d0, -0.00188d0/)
    real*8, parameter :: del_h = 0.074d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 450d0
    real*8, parameter :: a_grow = 0.05d0

    ! size of the human capital grid
    real*8 :: h_l, h_u
    real*8, parameter :: h_grow = 0.03d0

    ! maximum values of a and h by year
    real*8 :: a_max(JJ), h_max(JJ)

    ! net prices
    real*8 :: r, w_m, w_f

    ! transfer payments (old-age) and survival probabilities
    real*8 :: pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), ym_coh(JJ), yf_coh(JJ), l_coh(JJ)
    real*8 :: a_coh(JJ), h_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ)

    ! the permanent shock process
    real*8 :: dist_theta(NP), theta(NP)

    ! variables related to children
    integer :: nchild(JJ)
    real*8  :: pchild(JJ)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 3

    ! demographic and other model parameters
    real*8 :: eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), h(0:NH), aplus(JJ, 0:NA, 0:NH, NP, NS, NS)
    real*8 :: c(JJ, 0:NA, 0:NH, NP, NS, NS), l(JJ, 0:NA, 0:NH, NP, NS, NS)
    real*8 :: phi(JJ, 0:NA, 0:NH, NP, NS, NS), V(JJ, 0:NA, 0:NH, NP, NS, NS) = 0d0

    ! numerical variables
    real*8 :: RHS(0:NA, 0:NH, NP, NS, NS, 0:1), EV(0:NA, 0:NH, NP, NS, NS, 0:1)
    integer :: ij_com, ia_com, ih_com, ip_com, ism_com, isf_com, il_com
    real*8 :: cons_com

contains


    ! the first order condition
    function foc(x)

        implicit none
        real*8, intent(in) :: x
        real*8 :: foc, a_plus, varphi, tomorrow, wage_m, wage_f
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus  = x

        ! calculate wages (incl. child care costs) of men and women
        wage_m = w_m*eff(ij_com)*theta(ip_com)*eta(ism_com)
        wage_f = w_f*dble(il_com)*(exp(h(ih_com))*theta(ip_com)*eta(isf_com) - pchild(ij_com))

        ! calculate consumption
        cons_com = (1d0+r)*a(ia_com) + pen(ij_com) + wage_f + wage_m - a_plus

        ! calculate linear interpolation for future part of first order condition
        a_plus = max(a_plus, a_l)
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ial, ih_com, ip_com, ism_com, isf_com, il_com) + &
                  (1d0-varphi)*RHS(iar, ih_com, ip_com, ism_com, isf_com, il_com)

        ! calculate first order condition for consumption
        foc = cons_com/sqrt(2d0+dble(nchild(ij_com)))**(1d0-gamma) - tomorrow

    end function


    ! calculates marginal utility of consumption
    function margu(cons, ij)

        implicit none
        real*8, intent(in) :: cons
        integer, intent(in) :: ij
        real*8 :: margu

        margu = cons**(-1d0/gamma)/sqrt(2d0+dble(nchild(ij)))**egam

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, lpart, ij, ih, ip, ism, isf)

        implicit none
        integer, intent(in) :: lpart, ij, ih, ip, ism, isf
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
            valuefunc = max(varphi*EV(ial, ih, ip, ism, isf, lpart) + &
                (1d0-varphi)*EV(iar, ih, ip, ism, isf, lpart), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = (c_help/sqrt(2d0+dble(nchild(ij))))**egam/egam &
                        - nu*dble(lpart) + beta*psi(ij+1)*valuefunc

    end function

end module
