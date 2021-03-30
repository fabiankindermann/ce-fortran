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
    real*8, parameter :: nu = 0.36d0
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.40d0
    real*8, parameter :: delta = 0.019d0
    real*8, parameter :: rho = 0.95d0
    real*8, parameter :: sigma_eps = 0.000049d0

    ! government variables
    real*8, parameter :: gy = 0.15d0
    logical :: income_tax

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 40d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! time path of consumption and capital
    integer, parameter :: TT = 1000000
    real*8 :: k0 = 15d0
    real*8 :: c_t(0:TT), l_t(0:TT), k_t(0:TT), eta_t(0:TT), tauw_t(0:TT), tauy_t(0:TT)
    real*8 :: r_t(0:TT), w_t(0:TT), i_t(0:TT), g_t(0:TT), y_t(0:TT)
    integer :: is_t(0:TT), is_sim(0:TT)

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS), l(0:NK, NS), tau_w(0:NK, NS), tau_y(0:NK, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK, NS), coeff_c(NK+3, NS), coeff_l(NK+3, NS)
    real*8 :: RHS(0:NK, NS), coeff_r(NK+3, NS)
    real*8 :: con_lev, x_in
    logical :: check

    ! business cycle statistics
    real*8 :: E_c, E_g, E_i, E_k, E_l, E_r, E_w, E_y, E_tauw, E_tauy
    real*8 :: CV_c, CV_g, CV_i, CV_k, CV_l, CV_r, CV_w, CV_y, CV_tauw, CV_tauy
    real*8 :: Cor_c, Cor_g, Cor_i, Cor_k, Cor_l, Cor_r, Cor_w, Cor_tauw, Cor_tauy

    ! storage variables for welfare analysis
    real*8 :: WE_c(2), WE_leis(2), WCV_c(2), WCV_leis(2), WCor_leis(2)

    ! variables to communicate with function
    real*8 :: k_com, tau_w_com, tau_y_com
    integer :: is_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, gov, kplus, lab, c_act, r_act, w_act

        lab = max(min(x_in, 1d0-1d-10), 1d-10)

        ! calculate government expenditure
        gov = gy*exp(eta(is_com))*k_com**alpha*lab**(1d0-alpha)

        ! calculate factor prices
        r_act = exp(eta(is_com))*alpha*(k_com/lab)**(alpha-1d0)-delta
        w_act = exp(eta(is_com))*(1d0 - alpha)*(k_com/lab)**alpha

        ! calculate labor earnings or income tax
        if(income_tax)then
            tau_y_com = gov/(w_act*lab + r_act*k_com)
            tau_y_com = min(max(tau_y_com, -1d0), 1d0)
            tau_w_com = 0d0
        else
            tau_w_com = gov/(w_act*lab)
            tau_w_com = min(max(tau_w_com, -1d0), 1d0)
            tau_y_com = 0d0
        endif

        ! get actual consumption level
        c_act = nu/(1d0-nu)*(1d0-lab)*(1d0-tau_w_com-tau_y_com)*w_act

        ! calculate future capital
        kplus = (1d0-tau_w_com-tau_y_com)*w_act*lab + (1d0+r_act*(1d0-tau_y_com))*k_com - c_act

        ! evaluate spline function
        foc = spline_eval(max(min(kplus, k_u), k_l), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc = foc + (RHS(NK,is_com)-RHS(NK-1,is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif
        if(kplus < k_l)then
            foc = foc + (RHS(1,is_com)-RHS(0,is_com))/(k(1)-k(0))*(kplus-k_l)
        endif
        foc = foc**(-1d0/gamma)

        ! evaluate first order condition
        foc = nu*(c_act**nu*(1d0-lab)**(1d0-nu))**egam/c_act - foc

    end function

end module
