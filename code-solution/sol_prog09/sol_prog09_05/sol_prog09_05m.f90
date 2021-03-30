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

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.4d0
    real*8, parameter :: delta = 0.019d0
    real*8, parameter :: nu = 0.36d0

    ! government parameters
    real*8, parameter :: gy = 0.1d0

    ! numerical parameters
    real*8, parameter :: k_l = 10d0
    real*8, parameter :: k_u = 50d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 3000

    ! counter variables
    integer :: it, ik, iter

    ! variables for the time path
    integer, parameter :: TT = 360
    real*8 :: k0
    real*8 :: k_t(0:TT), l_t(0:TT), y_t(0:TT), c_t(0:TT), tau_r_t(0:TT), tau_w_t(0:TT)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK), l(0:NK), tau_r(0:NK), tau_w(0:NK), r(0:NK), w(0:NK)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK), coeff_c(NK+3), coeff_l(NK+3)
    real*8 :: coeff_tr(NK+3), coeff_tw(NK+3), coeff_r(NK+3), RHS(0:NK)
    real*8 :: con_lev, x_in
    logical :: check

    ! government and taxes
    real*8 :: gbar, tau_r_init, tau_w_init

    ! welfare levels
    real*8 :: V_stat, V_trans

    ! variables to communicate with function
    integer :: tax, ik_com
    real*8 :: k_com, tau_r_com, tau_w_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: r, w, foc, kplus, c_curr, lab

        ! labour supply
        lab = max(min(x_in, 1d0-1d-10), 1d-10)

        ! interest rate and wage income
        r = alpha*(lab/k_com)**(1d0-alpha) - delta
        w = (1d0-alpha)*(k_com/lab)**alpha

        ! budget balancing tax rate
        if(tax == 1)then
            tau_w_com = min((gbar - tau_r_init*r*k_com)/(w*lab), 1d0-1d-10)
        else
            tau_r_com = min((gbar - tau_w_init*w*lab)/(r*k_com), 1d0-1d-10)
        endif

        ! calculate optimal current consumption level
        c_curr = max(nu/(1d0-nu)*(1d0-tau_w_com)*w*(1d0-lab), 1d-10)

        ! derive next period capital level
        kplus = (1d0-delta)*k_com + k_com**alpha*lab**(1d0-alpha) - c_curr - gbar

        ! calculate future utility
        foc = spline_eval(max(min(kplus, k_u), k_l), coeff_r(:), k_l, k_u)
        if(kplus > k_u)then
            foc = foc + (RHS(NK)-RHS(NK-1))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif
        if(kplus < k_l)then
            foc = foc + (RHS(1)-RHS(0))/(k(1)-k(0))*(kplus-k_l)
        endif

        foc = max(foc, 1d-10)**(-1d0/gamma)

        ! get first-order condition
        foc = nu*(c_curr**nu*(1d0-lab)**(1d0-nu))**egam/c_curr - foc

    end function


    ! function that contains steady state conditions
    function foc_ss(x_in)

        implicit none
        real*8, intent(in)  :: x_in(:)
        real*8 :: x_in_temp(size(x_in, 1)), foc_ss(size(x_in, 1))

        ! restrict input variables
        x_in_temp(1) = x_in(1)
        x_in_temp(2) = max(min(x_in(2), 1d0-1d-10), 1d-10)
        x_in_temp(3) = max(min(x_in(3), 1d0-1d-10), -1d0+1d-10)

        ! solve first-order conditions
        if(tax == 1)then
            foc_ss(1) = nu/(1d0-nu)*(1d0-x_in_temp(3))*(1d0-alpha) &
                        *(x_in_temp(1)/x_in_temp(2))**alpha*(1d0-x_in_temp(2)) &
                        - x_in_temp(1)**alpha*x_in_temp(2)**(1d0-alpha) &
                        + delta*x_in_temp(1) + gbar
            foc_ss(2) = 1d0-beta*(1d0+(1d0-tau_r_init)*(alpha*(x_in_temp(1) &
                        /x_in_temp(2))**(alpha-1d0)-delta))
            foc_ss(3) = gbar - x_in_temp(3)*(1d0-alpha)*(x_in_temp(1) &
                        /x_in_temp(2))**alpha*x_in_temp(2) &
                        - tau_r_init*(alpha*(x_in_temp(1)/x_in_temp(2)) &
                        **(alpha-1d0)-delta)*x_in_temp(1)
        else
            foc_ss(1) = nu/(1d0-nu)*(1d0-tau_w_init)*(1d0-alpha) &
                        *(x_in_temp(1)/x_in_temp(2))**alpha*(1d0-x_in_temp(2)) &
                        - x_in_temp(1)**alpha*x_in_temp(2) &
                        **(1d0-alpha)+delta*x_in_temp(1)+gbar
            foc_ss(2) = 1d0-beta*(1d0+(1d0-x_in_temp(3))*(alpha*(x_in_temp(1) &
                        /x_in_temp(2))**(alpha-1d0)-delta))
            foc_ss(3) = gbar-tau_w_init*(1d0-alpha)*(x_in_temp(1) &
                        /x_in_temp(2))**alpha*x_in_temp(2) &
                        - x_in_temp(3)*(alpha*(x_in_temp(1)/x_in_temp(2)) &
                        **(alpha-1d0)-delta)*x_in_temp(1)
        endif

    end function

end module
