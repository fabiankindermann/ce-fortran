!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
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
    integer, parameter :: NA = 80

    ! number of persistent shock process values
    integer, parameter :: NP = 2

    ! number of transitory shock process values
    integer, parameter :: NS = 5

    ! number of points on the human capital grid
    integer, parameter :: NH = 60

    ! household preference parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: nu = 0.12d0
    real*8, parameter :: beta = 0.98d0

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

    ! the permanent shock process
    real*8 :: dist_theta(NP), theta(NP)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 3

    ! variables related to children
    integer :: nchild(JJ)
    real*8  :: pchild(JJ)

    ! net prices
    real*8 :: r, w_m, w_f

    ! transfer payments (old-age), survival probabilities, productivity
    real*8 :: pen(JJ), psi(JJ+1), eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), h(0:NH), aplus(JJ, 0:NA, 0:NH, NP, NS, NS)
    real*8 :: c(JJ, 0:NA, 0:NH, NP, NS, NS), l(JJ, 0:NA, 0:NH, NP, NS, NS)
    real*8 :: phi(JJ, 0:NA, 0:NH, NP, NS, NS), V(JJ, 0:NA, 0:NH, NP, NS, NS) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), ym_coh(JJ), yf_coh(JJ), l_coh(JJ)
    real*8 :: a_coh(JJ), h_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ)

    ! numerical variables
    real*8 :: EV(0:NA, 0:NH, NP, NS, NS, 0:1)
    integer :: ij_com, ia_com, ih_com, ip_com, ism_com, isf_com, il_com
    real*8 :: cons_com

contains


    ! calculates the value function
    function valuefunc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: valuefunc, a_plus, cons_lim, varphi, wage_m, wage_f, available, vplus
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate wages (incl. child care costs) of men and women
        wage_m = w_m*eff(ij_com)*theta(ip_com)*eta(ism_com)
        wage_f = w_f*dble(il_com)*(exp(h(ih_com))*theta(ip_com)*eta(isf_com) - pchild(ij_com))

        ! calculate available resources
        available = (1d0+r)*a(ia_com) + pen(ij_com) + wage_f + wage_m

        ! calculate consumption
        cons_com = available - a_plus
        cons_lim = max(cons_com, 1d-10)

        ! get interpolation weights
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate future utility
        vplus = (varphi*EV(ial, ih_com, ip_com, ism_com, isf_com, il_com) + &
          (1d0-varphi)*EV(iar, ih_com, ip_com, ism_com, isf_com, il_com))**egam/egam

        ! get value function
        valuefunc = -((cons_lim/sqrt(2d0+dble(nchild(ij_com))))**egam/egam &
                    - nu*dble(il_com) + beta*psi(ij_com+1)*vplus)

        ! penalty for negative consumption
        valuefunc = valuefunc + 100d0*abs(cons_com-cons_lim)

    end function

end module
