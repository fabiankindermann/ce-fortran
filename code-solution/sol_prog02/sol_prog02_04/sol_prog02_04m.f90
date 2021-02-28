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

    ! set parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: r = 0d0
    real*8, parameter :: beta = 1d0
    real*8, parameter :: w = 1d0

contains

    function utility(x_in)

        implicit none

        ! declaration of variables
        real*8, intent(in) :: x_in(:)
        real*8 :: utility
        real*8 :: c(3)

        ! define consumption levels
        c(2) = x_in(1)
        c(3) = x_in(2)
        c(1) = w + w/(1d0+r) - c(2)/(1d0+r) - c(3)/(1d0+r)**2d0

        ! set up utility we want to maximize
        utility = -(c(1)**egam + beta*c(2)**egam + beta**2d0*c(3)**egam)/egam

    end function

end module
