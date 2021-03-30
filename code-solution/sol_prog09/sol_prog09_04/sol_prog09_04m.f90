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
    real*8, parameter :: k0 = 10d0

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 100d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 3000

    ! counter variables
    integer :: it, ik, iter, count

    ! time path of consumption and capital
    integer, parameter :: TT = 360
    real*8 :: c_t(0:TT), k_t(0:TT), l_t(0:TT), y_t(0:TT)
    real*8 :: c_stat, k_stat, l_stat, y_stat

    ! value and policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK), l(0:NK), V(0:NK), V_new(0:NK)

    ! variables to numerically determine policy function
    real*8 :: coeff_c(NK+3), coeff_l(NK+3), coeff_V(NK+3)
    real*8 :: con_lev, x_in, fret
    logical :: check

    ! variables to communicate with function
    real*8 :: k_com

contains


    ! the functions that should be minimized
    function utility(x_in)

        use toolbox

        ! input variable labour supply
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, cons_lim, vplus
        real*8 :: kplus, kplus_lim, lab

        ! labour supply must be between 0 and 1
        lab = max(min(x_in, 1d0-1d-10), 1d-10)

        ! calculate current consumption
        cons = nu/(1d0-nu)*(1d0-lab)*(1d0-alpha)*(k_com/lab)**alpha
        cons_lim = max(cons, 1d-10)

        ! calculate future capital
        kplus = (1d0-delta)*k_com + k_com**alpha*lab**(1d0-alpha) - cons
        kplus_lim = max(min(kplus, k_u), k_l)

        ! calculate future utility
        vplus = max(spline_eval(kplus_lim, coeff_V, k_l, k_u), 1d-10)**egam/egam

        ! get utility function
        utility = -((cons**nu*(1d0-lab)**(1d0 - nu))**egam/egam + beta*vplus)

        ! penalty for negtive consumption and kplus outside of [k_l, k_u]
        utility = utility + abs(cons-cons_lim) + abs(kplus-kplus_lim)

    end function


    ! function that defines steady state
    function focs(x_in)

        implicit none
        real*8, intent(in)  :: x_in(:)
        real*8 :: k_temp, l_temp, focs(size(x_in, 1))

        k_temp = x_in(1)
        l_temp = max(min(x_in(2), 1d0-1d-10), 1d-10)

        focs(1) = nu/(1d0-nu)*(1d0-l_temp)/k_temp*(1d0-alpha)-l_temp/k_temp &
                  + delta*(l_temp/k_temp)**alpha
        focs(2) = 1d0 - beta*(1d0+alpha*(k_temp/l_temp)**(alpha-1d0)-delta)

    end function

end module
