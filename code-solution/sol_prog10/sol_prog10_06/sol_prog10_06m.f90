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

    ! year the household retires
    integer, parameter :: JR = 45

    ! year of definite infertilility
    integer, parameter :: JF = 25

    ! number of points on the asset grid
    integer, parameter :: NA = 50

    ! number of fertility shocks
    integer, parameter :: NF = 1

    ! number of children per year
    integer, parameter :: NC = 1

    ! number of total possible children
    integer, parameter :: NN = 10

    ! number of transitory shock process values
    integer, parameter :: NS = 7

    ! household preference and risk parameters
    real*8, parameter :: gamma = 0.50d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: chi = 5d0
    real*8, parameter :: echi = 1d0 - 1d0/chi
    real*8, parameter :: beta = 0.98d0
    real*8, parameter :: alpha = 0.12d0

    ! household risk process
    real*8, parameter :: sigma_eps = 0.22d0
    real*8, parameter :: rho       = 0.95d0

    ! probability to stay at home for child
    real*8, parameter :: p = 0.9d0

    ! household economies to time input
    real*8, parameter :: nu = 0.5d0

    ! raw time input per child
    real*8, parameter :: lc = 0.3d0

    ! numerical parameters
    real*8, parameter :: a_l    = 0d0
    real*8, parameter :: a_u    = 600d0
    real*8, parameter :: a_grow = 0.05d0

    ! probabilities to turn infertile
    real*8 :: pf(1:JJ, 0:NF, 0:NF)

    ! the transitory shock process
    real*8 :: pi(NS, NS), eta(NS)
    integer :: is_initial = 4

    ! net prices, transfer payments (old-age)
    real*8 :: r, w

    ! transfer payments (old-age), survival probabilities, productivity
    real*8 ::  pen(JJ), psi(JJ+1), eff(JJ)

    ! individual variables
    real*8 :: a(0:NA), aplus(JJ, 0:NA, 0:NN, NS, 0:NF), c(JJ, 0:NA, 0:NN, NS, 0:NF)
    real*8 :: chplus(JJ, 0:NA, 0:NN, NS, 0:NF), phi(JJ, 0:NA, 0:NN, NS, 0:NF)
    real*8 :: VC(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0, VN(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0
    real*8 :: V(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), ch_coh(JJ), n_ch_coh(JJ), a_coh(JJ), y_coh(JJ), v_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), age_ch_coh(0:NN), sum_phi(0:NN)

    ! numerical variables
    real*8 :: EVC(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0, EVN(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0
    real*8 :: RHS(JJ, 0:NA, 0:NN, NS, 0:NF) = 0d0
    integer :: ij_com, ia_com, in_com, is_com, ift_com
    integer :: chplus_com
    real*8 :: cons_com

contains


    ! the first order condition
    function foc(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, a_plus, varphi, tomorrow, wage, available
        integer :: ial, iar, in_h

        ! calculate tomorrows assets
        a_plus = x_in

        ! calculate the wage rate
        wage = w*eff(ij_com)*eta(is_com)*(1d0-in_com**nu*lc)

        ! calculate available resources
        available = (1d0+r)*a(ia_com) + wage + pen(ij_com)

        ! calculate consumption
        cons_com = available - a_plus

        ! calculate linear interpolation for future part of foc
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! get number of children
        in_h = in_com + chplus_com
        if(in_h > NN)in_h = NN

        tomorrow = varphi*RHS(ij_com+1, ial, in_h, is_com, ift_com) + &
                            (1d0-varphi)*RHS(ij_com+1, iar, in_h, is_com, ift_com)

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


    ! calculates the consumption value function
    function valuefunc_c(a_plus, cons, k_plus, ij, in, is, ift)

        implicit none
        integer, intent(in) :: k_plus, ij, in, is, ift
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc_c, varphi, c_help
        integer :: ial, iar, in_h

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! calculate linear interpolation for future part the value function
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! get number of children
        in_h = in + k_plus
        if(in_h > NN)in_h = NN

        ! calculate tomorrow's part of the value function
        valuefunc_c = 0d0
        if(ij < JJ)then
            valuefunc_c = max(varphi*EVC(ij+1, ial, in_h, is, ift) + &
                         (1d0-varphi)*EVC(ij+1, iar, in_h, is, ift), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc_c = c_help**egam/egam + beta*psi(ij+1)*valuefunc_c

    end function

    ! calculates the children value function
    function valuefunc_n(a_plus, k_plus, ij, in, is, ift)

        implicit none
        integer, intent(in) :: k_plus, ij, in, is, ift
        real*8, intent(in) :: a_plus
        real*8 :: valuefunc_n, varphi
        integer :: ial, iar, in_h

        ! calculate linear interpolation for future part the value function
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! get number of children
        in_h = in + k_plus
        if(in_h > NN)in_h = NN

        ! calculate tomorrow's part of the value function
        valuefunc_n = 0d0
        if(ij < JJ)then
            valuefunc_n = max(varphi*EVN(ij+1, ial, in_h, is, ift) + &
                         (1d0-varphi)*EVN(ij+1, iar, in_h, is, ift), 1d-10)**echi/echi
        endif

        ! add todays part and discount
        valuefunc_n = alpha*in**echi/echi + beta*psi(ij+1)*valuefunc_n

    end function

end module
