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
    real*8, parameter  :: psi(2:3) = (/0.8d0, 0.5d0/)
    real*8, parameter  :: mu_w = 1d0
    real*8, parameter  :: sig_w = 0.9d0
    integer, parameter :: n_w = 5

    real*8, parameter  :: R = 1d0
    real*8, parameter  :: beta = 1d0
    real*8, parameter  :: gamma = 0.5d0
    real*8, parameter  :: egam = 1d0-1d0/gamma
    real*8, parameter  :: nu = 0.05d0
    real*8, parameter  :: pen = mu_w

    real*8 :: w(n_w), weight_w(n_w), p_a(2)
    real*8 :: a(4, n_w), omega(2, n_w), c(3, n_w)
    real*8 :: wag(3, n_w), inc(3, n_w), sav(3, n_w), alp(3, n_w), b(2:4, n_w)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility
        integer :: iw, ic

        ! savings
        a(1, :) = 0d0
        a(2, :) = x(1)
        omega(1, :) = x(2)
        ic = 3
        do iw = 1, n_w
            a(3, iw) = x(ic)
            omega(2, iw) = x(ic + 1)
            a(4, iw) = x(ic + 2)
            ic = ic + 3
        enddo

        ! consumption (insure consumption > 0)
        c(1, :) = mu_w - a(2, 1)
        c(2, :) = R*(1d0 - omega(1, 1))*a(2, 1) + omega(1, 1)*a(2, 1)/p_a(1) + &
                 w(:) - a(3,:)
        c(3,:) = R*(1d0-omega(2, :))*a(3, :) + omega(2, :)*a(3, :)/p_a(2) + &
                 omega(1, 1)*a(2, 1)/p_a(1) + pen - a(4, :)
        c = max(c, 1d-10)

        ! bequests (insure bequests > 0)
        b(2, :) = R*(1d0 - omega(1, 1))*a(2, 1)
        b(3, :) = R*(1d0 - omega(2, :))*a(3, :)
        b(4, :) = R*a(4, :)
        b = max(b, 1d-10)

        ! expected utility of period 3
        utility = 0d0
        do iw = 1, n_w
            utility = utility + weight_w(iw)*beta**2*psi(2)*psi(3)*         &
                                            (c(3, iw)**egam + beta*nu*b(4, iw)**egam)
        enddo

        ! expected utility of period 2
        do iw = 1, n_w
            utility = utility + weight_w(iw)*beta*psi(2)*       &
                               (c(2, iw)**egam + beta*(1d0-psi(3))*nu*b(3, iw)**egam)
        enddo

        ! add first period utility
        utility = -(c(1, 1)**egam + beta*(1d0-psi(2))*nu*b(2, 1)**egam + utility)/egam

    end function


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
