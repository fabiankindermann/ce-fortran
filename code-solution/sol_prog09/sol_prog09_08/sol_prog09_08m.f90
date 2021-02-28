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

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: nu = 0.36d0
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.40d0
    real*8, parameter :: delta = 0.019d0
    real*8, parameter :: rho = 0.95d0
    real*8, parameter :: sigma_eps = 0.000049d0

    ! numerical parameters
    real*8, parameter :: k_l = 10d0
    real*8, parameter :: k_u = 40d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, is, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 100000
    real*8 :: k0 = 15d0
    real*8 :: c_t(0:TT), l_t(0:TT), k_t(0:TT), eta_t(0:TT)
    real*8 :: r_t(0:TT), w_t(0:TT), i_t(0:TT), y_t(0:TT)
    integer :: is_t(0:TT)

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy and value functions
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS), l(0:NK, NS), V(0:NK, NS)

    ! variables to numerically determine value and policy function
    real*8 :: V_new(0:NK, NS), c_new(0:NK, NS), coeff_c(NK+3, NS), coeff_l(NK+3, NS)
    real*8 :: EV(0:NK, NS), coeff_ev(NK+3, NS)
    real*8 :: RHS(0:NK, NS), coeff_r(NK+3, NS)
    real*8 :: con_lev, x_in, y_in(2), fret
    logical :: check

    ! business cycle statistics
    real*8 :: E_c, E_i, E_k, E_l, E_r, E_w, E_y
    real*8 :: CV_c, CV_i, CV_k, CV_l, CV_r, CV_w, CV_y
    real*8 :: Cor_c, Cor_i, Cor_k, Cor_l, Cor_r, Cor_w

    ! variables to communicate with function
    real*8 :: k_com
    integer :: is_com

contains


    ! the value function
    function utility(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, evplus, kplus, kplus_lim, lab

        ! labor supply must be between 0 and 1
        lab = max(min(x_in, 1d0-1d-10), 1d-10)

        ! calculate current consumption
        cons = max((nu/(1d0-nu))*(1d0-lab)*exp(eta(is_com))*(1d0-alpha)* &
                (k_com/lab)**alpha, 1d-10)

        ! calculate future capital
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha* &
                    lab**(1d0-alpha) - cons
        kplus_lim = max(min(kplus,k_u),k_l)

        ! calculate future utility
        evplus = max(spline_eval(kplus_lim, coeff_ev(:,is_com), &
                    k_l, k_u), 1d-10)**egam/egam

        ! get utility function
        utility = -((cons**nu*(1d0-lab)**(1d0-nu))**egam/egam + beta*evplus)

        ! penalize for kplus values outside of (k_l, k_u)
        utility = utility + 100d0*abs(kplus-kplus_lim)

    end function


    ! the first order condition
    function foc(y_in)

        use toolbox

        implicit none
        real*8, intent(in) :: y_in(:)
        real*8 :: kplus, y_in_temp(size(y_in, 1))
        real*8 :: foc(size(y_in, 1))

        y_in_temp(1) = max(y_in(1), 1d-10)
        y_in_temp(2) = max(min(y_in(2), 1d0-1d-10), 1d-10)

        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha*y_in_temp(2)**(1d0-alpha) - y_in_temp(1)
        kplus = max(min(kplus, k_u), k_l)

        ! calculate optimal actual consumption level
        foc(1) = y_in_temp(1) - nu/(1d0-nu)*(1d0-y_in_temp(2))&
                    *exp(eta(is_com))*(1d0-alpha)*(k_com/y_in_temp(2))**alpha
        foc(2) = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        foc(2) = nu*(y_in_temp(1)**nu*(1d0-y_in_temp(2))**(1d0-nu))**egam/y_in_temp(1) - foc(2)**(-1d0/gamma)

    end function


    ! the first order condition in consumption values
    function foc2(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc2, kplus, c_act

        ! calculate optimal actual consumption level
        c_act = nu/(1d0-nu)*(1d0-x_in)*exp(eta(is_com))*(1d0-alpha)*(k_com/x_in)**alpha

        ! calculate future capital
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha*x_in**(1d0-alpha) - c_act

        ! evaluate spline function
        foc2 = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc2 = foc2 + (RHS(NK, is_com) - RHS(NK-1, is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif
        foc2 = foc2**(-1d0/gamma)

        ! evaluate first order condition
        foc2 = c_act - nu*(c_act**nu*(1d0-x_in)**(1d0-nu))**egam/foc2

    end function

end module
