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
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.40d0
    real*8, parameter :: k0 = 0.1d0

    ! numerical parameters
    real*8, parameter :: k_l = 0.1d0
    real*8, parameter :: k_u = 0.4d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 30
    real*8 :: c_t(0:TT), k_t(0:TT), y_t(0:TT)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK), V(0:NK)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK), V_new(0:NK), coeff_c(NK+3), coeff_V(NK+3)
    real*8 :: con_lev, x_in, fret
    logical :: check

    ! variables to communicate with function
    real*8 :: k_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, kplus, cplus

        ! future capital
        kplus = k_com**alpha - x_in

        ! calculate future consumption
        cplus = spline_eval(kplus, coeff_c, k_l, k_u)

        ! get first order condition
        foc = 1d0/x_in - beta*alpha*kplus**(alpha-1d0)/cplus

    end function

    ! the function that should be minimized
    function utility(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, vplus

        ! calculate consumption
        cons = k_com**alpha - x_in

        ! calculate future utility
        vplus = spline_eval(x_in, coeff_V, k_l, k_u)

        ! get utility function
        if(cons < 1d-10)then
            utility = 1d10*(1d0 + abs(cons))
        else
            utility = -(log(cons) + beta*vplus)
        endif

    end function

end module
