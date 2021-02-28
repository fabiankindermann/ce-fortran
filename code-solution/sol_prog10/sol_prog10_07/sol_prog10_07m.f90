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
    real*8, parameter :: nu = 0.335d0
    real*8, parameter :: beta = 0.98d0

    ! household risk process
    real*8, parameter :: sigma_theta = 0.242d0
    real*8, parameter :: sigma_eps   = 0.022d0
    real*8, parameter :: rho         = 0.985d0

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
    real*8 :: pen(JJ), psi(JJ+1), eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, NP, NS)
    real*8 :: c(JJ, 0:NA, NP, NS), l(JJ, 0:NA, NP, NS)
    real*8 :: phi(JJ, 0:NA, NP, NS), V(JJ, 0:NA, NP, NS) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), l_coh(JJ), h_coh(JJ), a_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_l(JJ), cv_h(JJ), corr_hl(JJ)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, NP, NS), EV(JJ, 0:NA, NP, NS)
    integer :: ij_com, ia_com, ip_com, is_com
    real*8 :: cons_com, lab_com

contains


    ! the value function
    function valuefunc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: valuefunc, a_plus, cons_lim, varphi, wage, available, vplus
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! calculate available resources
        available = (1d0+r)*a(ia_com) + pen(ij_com)

        ! determine labor
        if(ij_com < JR)then
            lab_com = min(max(nu + (1d0-nu)*(a_plus - available)/wage, 0d0), 1d0-1d-10)
        else
            lab_com = 0d0
        endif

        ! calculate consumption
        cons_com = available + wage*lab_com - a_plus
        cons_lim = max(cons_com, 1d-10)

        ! get interpolation weights
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate future utility
        vplus = (varphi*EV(ij_com+1, ial, ip_com, is_com) + &
                 (1d0-varphi)*EV(ij_com+1, iar, ip_com, is_com))**egam/egam

        ! get value function
        valuefunc = -((cons_com**nu*(1d0-lab_com)**(1d0-nu))**egam/egam + &
                      beta*psi(ij_com+1)*vplus)

        ! penalty for negative consumption
        valuefunc = valuefunc + 100d0*abs(cons_com - cons_lim)

    end function


    ! calculates marginal utility of consumption
    function margu(cons, lab)

        implicit none
        real*8, intent(in) :: cons, lab
        real*8 :: margu

        margu = nu*(cons**nu*(1d0-lab)**(1d0-nu))**egam/cons

    end function


    ! calculates Euler equation error
    function eer(a_err, a_plus)

        implicit none
        real*8, intent(in) :: a_err, a_plus
        real*8 :: eer, cons, lab, available, varphi, tomorrow, wage
        integer :: ial, iar

        ! calculate the wage rate
        wage = w*eff(ij_com)*theta(ip_com)*eta(is_com)

        ! calculate available resources
        available = (1d0+r)*a_err + pen(ij_com)

        ! determine labor
        if(ij_com < JR)then
            lab = min(max(nu + (1d0-nu)*(a_plus-available)/wage, 0d0), 1d0-1d-10)
        else
            lab = 0d0
        endif

        ! calculate consumption
        cons = available + wage*lab - a_plus
        cons = max(cons, 1d-10)

        ! get interpolation weights
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)
        tomorrow = varphi*RHS(ij_com+1, ial, ip_com, is_com) + &
                   (1d0-varphi)*RHS(ij_com+1, iar, ip_com, is_com)
        tomorrow = tomorrow**(-1d0/gamma)

        ! calculate Euler equation residual (eer)
        eer = 1d0 - margu(cons, lab)/tomorrow

    end function

end module
