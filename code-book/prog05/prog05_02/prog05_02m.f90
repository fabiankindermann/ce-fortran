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
    real*8, parameter  :: mu_w = 1d0
    real*8, parameter  :: sig_w = 0.0d0
    integer, parameter :: n_w = 5

    real*8, parameter  :: R = 1d0
    real*8, parameter  :: beta = 1d0
    real*8, parameter  :: gamma = 0.5d0
    real*8, parameter  :: egam = 1d0-1d0/gamma

    real*8 :: w(n_w), weight_w(n_w)
    real*8 :: a(3,n_w), c(3,n_w)
    real*8 :: wag(3,n_w), inc(3,n_w), sav(3,n_w)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility
        integer :: iw

        ! savings
        a(1, :) = 0d0
        a(2, :) = x(1)
        a(3, :) = x(2:1+n_w)

        ! consumption (ensure consumption > 0)
        c(1,:) = mu_w - a(2,1)
        c(2,:) = R*a(2,1) + w(:) - a(3,:)
        c(3,:) = R*a(3,:)
        c = max(c, 1d-10)

        ! expected utility of periods 2 and 3
        utility = 0d0
        do iw = 1, n_w
            utility = utility+weight_w(iw)*(c(2,iw)**egam + beta*c(3,iw)**egam)
        enddo
        utility = -(c(1,1)**egam + beta*utility)/egam

    end function


    ! calculates the expected value of x
    function E(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: E
        integer :: iw

        ! calculate expected value
        E = 0d0
        do iw = 1, n_w
            E = e + x(iw)*weight_w(iw)
        enddo

    end function


    ! calculates the standard deviation of x
    function Std(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: Std
        real*8 :: E
        integer :: iw

        ! calculate expected value
        E = 0d0
        do iw = 1, n_w
            E = E + x(iw)*weight_w(iw)
        enddo

        ! calculate standard deviation
        Std = 0d0
        do iw = 1, n_w
            Std = Std + x(iw)**2*weight_w(iw)
        enddo
        Std = sqrt(max(Std-E**2, 0d0))

    end function

end module
