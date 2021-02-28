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
    real*8, parameter :: w = 1d0
    real*8, parameter :: R = 1.0d0
    real*8, parameter :: beta = 1d0
    real*8, parameter :: theta = 0.5d0
    real*8, parameter :: nu = 0.5d0
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: delh = 1.0d0
    real*8, parameter :: a_low = 0.5d0

    real*8 :: a(3), c(3), ah, p, u(3)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility

        ! savings
        a(1) = 0d0
        ah = x(1)
        a(2) = x(2)
        a(3) = x(3)

        p = R - 1d0 + delh

        ! consumption (insure consumption > 0)
        c(1) = w - a(2) - ah
        c(2) = R*a(2) + w - a(3) - delh*ah
        c(3) = R*a(3) + (1d0 - p - delh)*ah
        c = max(c, 1d-10)
        ah = max(ah, 1d-10)

        ! utility function
        u(1) = (theta*c(1)**nu+(1d0-theta)*ah**nu)**(egam/nu)/egam
        u(2) = (theta*c(2)**nu+(1d0-theta)*ah**nu)**(egam/nu)/egam
        u(3) = (theta*c(3)**nu+(1d0-theta)*ah**nu)**(egam/nu)/egam

        utility = -(u(1) + beta*u(2) + beta**2*u(3))

    end function

end module
