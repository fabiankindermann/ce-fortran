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
    logical :: variable_g
    real*8 :: g_level

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 40d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! time path of consumption and capital
    integer, parameter :: TT = 1000000
    real*8 :: k0 = 15d0
    real*8 :: c_t(0:TT), l_t(0:TT), k_t(0:TT), eta_t(0:TT)
    real*8 :: r_t(0:TT), w_t(0:TT), i_t(0:TT), g_t(0:TT), y_t(0:TT)
    integer :: is_t(0:TT), is_sim(0:TT)

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS), l(0:NK, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK, NS), coeff_c(NK+3, NS), coeff_l(NK+3, NS)
    real*8 :: RHS(0:NK, NS), coeff_r(NK+3, NS)
    real*8 :: con_lev, x_in
    logical :: check

    ! business cycle statistics
    real*8 :: E_c, E_g, E_i, E_k, E_l, E_r, E_w, E_y
    real*8 :: CV_c, CV_g, CV_i, CV_k, CV_l, CV_r, CV_w, CV_y
    real*8 :: Cor_c, Cor_g, Cor_i, Cor_k, Cor_l, Cor_r, Cor_w

    ! storage variables for welfare analysis
    real*8 :: WE_c(2), WE_leis(2), WCV_c(2), WCV_leis(2), WCor_leis(2)

    ! variables to communicate with function
    real*8 :: k_com
    integer :: is_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, kplus, c_act, gov

        ! calculate optimal actual consumption level
        c_act = nu/(1d0-nu)*(1d0-x_in)*exp(eta(is_com))*(1d0-alpha)*(k_com/x_in)**alpha

        ! calculate government expenditure
        if(variable_g)then
            gov = gy*exp(eta(is_com))*k_com**alpha*x_in**(1d0-alpha)
        else
            gov = g_level
        endif

        ! calculate future capital
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha*x_in**(1d0-alpha) - c_act - gov

        ! evaluate spline function
        foc = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc = foc + (RHS(NK,is_com)-RHS(NK-1,is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif
        foc = foc**(-1d0/gamma)

        ! evaluate first order condition
        foc = nu*(c_act**nu*(1d0-x_in)**(1d0-nu))**egam/c_act - foc

    end function


    ! the first order condition in consumption values
    function foc2(x_in)

        ! modules to be used
        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc2, kplus, c_act, gov

        ! calculate optimal actual consumption level
        c_act = nu/(1d0-nu)*(1d0-x_in)*exp(eta(is_com))*(1d0-alpha)*(k_com/x_in)**alpha

        ! calculate government expenditure
        if(variable_g)then
            gov = gy*exp(eta(is_com))*k_com**alpha*x_in**(1d0-alpha)
        else
            gov = g_level
        endif

        ! calculate future capital
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha*x_in**(1d0-alpha) - c_act - gov

        ! evaluate spline function
        foc2 = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc2 = foc2 + (RHS(NK,is_com)-RHS(NK-1,is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif
        foc2 = foc2**(-1d0/gamma)

        ! evaluate first order condition
        foc2 = c_act - nu*(c_act**nu*(1d0-x_in)**(1d0-nu))**egam/foc2

    end function

end module
