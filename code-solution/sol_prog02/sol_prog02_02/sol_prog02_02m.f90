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

    ! declaration of variables
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 1d0
    real*8, parameter :: r = 0d0
    real*8, parameter :: w = 1d0

contains

    ! the first order conditions from the lagrangian
    function foc(x_in)

        implicit none

        real*8, intent(in) :: x_in(:)
        real*8 :: foc(size(x_in, 1))

        ! set up equation system to solve
        foc(1) = x_in(1)**(-1d0/gamma) - x_in(3)
        foc(2) = beta*x_in(2)**(-1d0/gamma) - x_in(3)/(1d0+r)
        foc(3) = w - x_in(1) - x_in(2)/(1d0+r)

    end function

    ! the utility function that should be maximized
    function utility(x_in)

        implicit none

        real*8, intent(in) :: x_in
        real*8 :: utility

        ! set up utility we want to maximize
        utility = -(x_in**egam/egam + beta*((w-x_in)*(1d0+r))**egam/egam)

    end function

end module
