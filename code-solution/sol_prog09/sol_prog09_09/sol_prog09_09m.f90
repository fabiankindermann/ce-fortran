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

    ! model parameters as in Hansen(1985)
    real*8, parameter :: A = 2d0
    real*8, parameter :: h0 = 0.53d0
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.36d0
    real*8, parameter :: delta = 0.025d0
    real*8, parameter :: rho = 0.95d0
    real*8, parameter :: sigma_eps = (0.712d-2)**2d0

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 250d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, is, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 100000
    real*8 :: k0 = 100d0
    real*8 :: c_t(0:TT), kappa_t(0:TT), k_t(0:TT), eta_t(0:TT)
    real*8 :: r_t(0:TT), w_t(0:TT), i_t(0:TT), y_t(0:TT)
    integer :: is_t(0:TT)

    ! Analytical Steady State variables
    real*8 :: kbar, hbar, ybar, cbar, rbar, wbar, ibar

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS), kappa(0:NK, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK, NS), coeff_c(NK+3, NS), coeff_kappa(NK+3, NS)
    real*8 :: RHS(0:NK, NS), coeff_r(NK+3, NS)
    real*8 :: con_lev, x_in
    logical :: check

    ! business cycle statistics
    real*8 :: E_c, E_i, E_kappa, E_k, E_r, E_w, E_y
    real*8 :: CV_c, CV_i, CV_k, CV_kappa, CV_r, CV_w, CV_y
    real*8 :: Cor_c, Cor_i, Cor_k, Cor_kappa, Cor_r, Cor_w

    ! variables to communicate with function
    real*8 :: k_com
    integer :: is_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, c_act, kappa_act, kplus

        ! restrict kappa to be in [0,1]
        kappa_act = min(max(x_in, 1d-10), 1d0-1d-10)

        ! calculate optimal actual consumption level
        c_act = - (1d0-alpha)*exp(eta(is_com))*(k_com/(kappa_act*h0))**alpha*h0/(A*log(1d0-h0))

        ! calculate future capital
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha*(kappa_act*h0)**(1d0-alpha) - c_act

        ! evaluate spline function
        foc = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc = foc + (RHS(NK, is_com)-RHS(NK-1, is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif

        ! evaluate first order condition
        foc = 1d0/c_act - foc

    end function

end module
