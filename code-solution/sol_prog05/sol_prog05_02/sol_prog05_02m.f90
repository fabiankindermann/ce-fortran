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
    real*8, parameter :: mu = 0.5d0
    real*8, parameter :: sigma = 0.5d0
    real*8, parameter :: w = 1d0
    real*8, parameter :: R = 1d0
    real*8, parameter :: beta = 1d0
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    integer, parameter :: NH = 20

    real*8 :: a(3), c(2), c3(NH)
    real*8 :: weight_h(NH), hc(NH)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility
        integer ::   ih

        ! savings
        a(1) = 0d0
        a(2) = x(1)
        a(3) = x(2)

        ! consumption (insure consumption > 0)
        c(1) = w - a(2)
        c(2) = R*a(2) + w - a(3)
        c3 = R*a(3) - hc
        c = max(c, 1d-10)
        c3 = max(c3, 1d-10)

        ! expected utility of period 3
        utility = 0d0
        do ih = 1, NH
            utility = utility + weight_h(ih)*c3(ih)**egam
        enddo

        ! utility function
        utility = -(c(1)**egam + beta*c(2)**egam + beta**2*utility)/egam

    end function

    ! calculates the expected value of x
    function E(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: E
        integer :: ih

        ! calculate expected value
        E = 0d0
        do ih = 1, NH
            E = E + x(ih)*weight_h(ih)
        enddo

    end function


    ! calculates the standard deviation of x
    function Std(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: Std
        real*8 :: E
        integer :: ih

        ! calculate expected value
        E = 0d0
        do ih = 1, NH
            E = E + x(ih)*weight_h(ih)
        enddo

        ! calculate standard deviation
        Std = 0d0
        do ih = 1, NH
            Std = Std + x(ih)**2*weight_h(ih)
        enddo
        Std = sqrt(max(Std-E**2, 0d0))

    end function

end module
