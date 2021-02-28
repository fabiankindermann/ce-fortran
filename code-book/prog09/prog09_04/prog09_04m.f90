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
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.40d0
    real*8, parameter :: delta = 0.019d0
    real*8, parameter :: k0 = 79d0
    real*8, parameter :: rho = 0.95d0
    real*8, parameter :: sigma_eps = 0.000049d0

    ! numerical parameters
    real*8, parameter :: k_l = 50d0
    real*8, parameter :: k_u = 110d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, is, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 5000
    real*8 :: c_t(0:TT), k_t(0:TT), eta_t(0:TT)
    integer :: is_t(0:TT)

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK, NS), coeff_c(NK+3, NS)
    real*8 :: RHS(0:NK, NS), coeff_r(NK+3, NS)
    real*8 :: con_lev, x_in
    logical :: check

    ! variables to communicate with function
    real*8 :: k_com
    integer :: is_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, kplus

        ! future capital stock
        kplus = (1d0-delta)*k_com + exp(eta(is_com))*k_com**alpha - x_in

        ! calculate future expected marginal utility
        foc = spline_eval(min(kplus, k_u), coeff_r(:, is_com), k_l, k_u)
        if(kplus > k_u)then
            foc = foc + (RHS(NK,is_com)-RHS(NK-1,is_com))/(k(NK)-k(NK-1))*(kplus-k_u)
        endif

        ! get first order condition
        foc = x_in - foc

    end function

end module
