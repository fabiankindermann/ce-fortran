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
    real*8, parameter :: beta = 0.975d0
    real*8, parameter :: r = 0.02d0
    real*8, parameter :: w = 1d0
    real*8, parameter :: a0 = 0d0
    real*8, parameter :: a_borrow = 5d0

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! time path of consumption and resource
    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT), a_t(0:TT)

    ! value and policy function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA), V(0:NA)
    real*8 :: a_l, a_u

    ! variables to numerically determine value and policy function
    real*8 :: V_new(0:NA), coeff_V(NA+3), coeff_c(NA+3)
    real*8 :: con_lev, x_in, fret

    ! variables to communicate with function
    integer :: ia_com

contains


    ! the function that should be minimized
    function utility(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: utility, cons, vplus

        ! calculate consumption
        cons = (1d0+r)*a(ia_com) + w - x_in

        ! calculate future utility
        vplus = max(spline_eval(x_in, coeff_V, a_l, a_u), 1d-10)**egam/egam

        ! get utility function
        if(cons < 1d-10)then
            utility = -1d-10**egam/egam*(1d0+abs(cons))
        else
            utility = - (cons**egam/egam + beta*vplus)
        endif

    end function

end module
