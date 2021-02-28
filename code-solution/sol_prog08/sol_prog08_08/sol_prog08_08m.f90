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
    real*8 :: a(0:NA), c(0:NA)
    real*8 :: a_l, a_u

    ! variables to numerically determine value and policy function
    real*8 :: c_new(0:NA), coeff_c(NA+3)
    real*8 :: con_lev, x_in, fret
    logical :: check

    ! variables to communicate with function
    integer :: ia_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        ! variable declarations
        real*8, intent(in) :: x_in
        real*8 :: foc, cplus

        ! calculate future utility
        cplus = spline_eval(min(x_in, a_u), coeff_c, a_l, a_u)

        ! extrapolate if p goes beyond a_u
        if(x_in > a_u)then
            cplus = cplus + (x_in-a_u)*(c(NA)-c(NA-1))/(a(NA)-a(NA-1))
        endif

        ! get utility function
        foc = (1d0+r)*a(ia_com) + w - x_in - (beta*(1d0+r))**(-gamma)*cplus

    end function

end module
