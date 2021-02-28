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
    real*8, parameter :: w = 1d0
    real*8, parameter :: R = 1d0
    real*8, parameter :: beta = 1d0
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma

    real*8 :: a(3), c(3)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility

        ! savings
        a(1) = 0d0
        a(2:3) = x

        ! consumption (insure consumption > 0)
        c(1) = w - a(2)
        c(2) = R*a(2) + w - a(3)
        c(3) = R*a(3)
        c = max(c, 1d-10)

        ! utility function
        utility = -(c(1)**egam + beta*c(2)**egam + beta**2*c(3)**egam)/egam

    end function

end module
