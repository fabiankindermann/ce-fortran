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
    integer :: it, ik, is, is_p, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 5000
    real*8 :: c_t(0:TT), k_t(0:TT), eta_t(0:TT)
    integer :: is_t(0:TT)

    ! the shock process
    integer, parameter :: NS = 21
    real*8 :: pi(NS, NS), eta(NS)

    ! policy and value functions
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK, NS), X(0:NK, NS), V(0:NK, NS)

    ! variables to numerically determine policy and value function
    real*8 :: VX(0:NK, NS), VX_endog(0:NK, NS), dV(0:NK, NS), EV(0:NK, NS)
    real*8 :: V_new(0:NK, NS), coeff_c(NK + 3, NS), coeff_EV(NK + 3, NS)
    real*8 :: k_endog(0:NK, NS), c_endog(0:NK, NS), X_endog(0:NK, NS)
    real*8 :: con_lev, x_in, fret
    logical :: check

    ! variables to communicate with function
    real*8 :: k_com, X_com
    integer :: is_com

contains


    ! the functions that should be minimized
    function utility(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, evplus, kplus

        ! calculate tomorrows capital
        kplus = max(min(x_in, k_u), k_l)

        ! calculate consumption
        cons = max(X_com - kplus, 1d-10)

        ! calculate future utility
        evplus = max(spline_eval(kplus, coeff_EV(:, is_com), &
                                            k_l, k_u), 1d-10)**egam/egam
        ! get utility function
        utility = -(cons**egam/egam + beta*evplus)

    end function


    ! derive the endogenous capital stock
    function fkend(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: fkend

        ! get first order condition
        fkend = X_com - exp(eta(is_com))*x_in**alpha - (1d0-delta)*x_in

    end function


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, kplus, cplus
        integer :: is_p

        ! future capital stock
        kplus = X_com - x_in

        ! calculate future expected marginal utility
        foc = 0d0
        do is_p = 1, NS
            cplus = spline_eval(min(kplus, k_u), coeff_c(:, is_p), k_l, k_u)
            if(kplus > k_u)then
                cplus = cplus + (c(NK,is_p) - c(NK-1,is_p))/(k(NK)-k(NK-1))*(kplus-k_u)
            endif
            foc = foc + pi(is_com, is_p)* &
                     (1d0+exp(eta(is_p))*alpha*kplus**(alpha-1d0)-delta)*cplus**(-1d0/gamma)
        enddo

        ! get first order condition
        foc = x_in - (beta*foc)**(-gamma)

    end function

end module
