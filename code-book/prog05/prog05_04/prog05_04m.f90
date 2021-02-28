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

    real*8, parameter  :: mu_R = 1.22d0
    real*8, parameter  :: sig_R = 0.0d0
    integer, parameter :: n_R = 5
    real*8, parameter  :: rho_wR = 0d0

    real*8, parameter  :: Rf = 1d0
    real*8, parameter  :: beta = 1d0
    real*8, parameter  :: gamma = 0.5d0
    real*8, parameter  :: egam = 1d0 - 1d0/gamma

    real*8 :: wR(n_w*n_R, 2), weight_wR(n_w*n_R)
    real*8 :: R(n_R), weight_R(n_R), hga
    real*8 :: a(3,n_w*n_R), omega(2,n_w*n_R), c(3,n_w*n_R,n_R)

contains


    ! utility function of the household
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility
        real*8 :: prob
        integer :: iw, ir2, ir3, ic, iwR

        ! savings
        a(1,:) = 0d0
        a(2,:) = x(1)
        omega(1,:) = x(2)
        ic = 3
        iwR = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                a(3, iwR) = x(ic)
                omega(2, iwR) = x(ic+1)
                ic = ic+2
                iwR = iwR+1
            enddo
        enddo

        ! consumption
        c(1,:,:) = mu_w - a(2,1)
        iwR = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                c(2,iwR,:) = (Rf+omega(1,1)*(wR(iwR,2)-Rf))*a(2,1) + &
                              wR(iwR,1) - a(3,iwR)
                do ir3 = 1, n_R
                    c(3,iwR,ir3) = (Rf+omega(2,iwR)*(R(ir3)-Rf))*a(3,iwR)
                enddo
                iwR = iwR+1
            enddo
        enddo
        c = max(c, 1d-10)

        ! expected utility of periods 2 and 3
        utility = 0d0
        iwR = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                do ir3 = 1, n_R
                    prob = weight_wR(iwR)*weight_R(ir3)
                    utility = utility+prob*(c(2,iwR,1)**egam+beta*c(3,iwR,ir3)**egam)
                enddo
                iwR = iwR+1
            enddo
        enddo
        utility = -(c(1,1,1)**egam + beta*utility)/egam

    end function


    ! calculates the expected value of x
    function E(x)

        implicit none
        real*8, intent(in) :: x(:,:)
        real*8 :: E
        integer :: iw, ir2, ir3, ic

        ! calculate expected value
        E = 0d0
        ic = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                do ir3 = 1, n_R
                    E = E + x(ic, ir3)*weight_wR(ic)*weight_R(ir3)
                enddo
                ic = ic+1
            enddo
        enddo

    end function


    ! calculates the standard deviation of x
    function Std(x)

        implicit none
        real*8, intent(in) :: x(:,:)
        real*8 :: Std
        real*8 :: E
        integer :: iw, ir2, ir3, ic

        ! calculate expected value
        E = 0d0
        ic = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                do ir3 = 1, n_R
                    E = E + x(ic, ir3)*weight_wR(ic)*weight_R(ir3)
                enddo
                ic = ic+1
            enddo
        enddo

        ! calculate standard deviation
        Std = 0d0
        ic = 1
        do iw = 1, n_w
            do ir2 = 1, n_R
                do ir3 = 1, n_R
                    Std = Std + x(ic, ir3)**2*weight_wR(ic)*weight_R(ir3)
                enddo
                ic = ic+1
            enddo
        enddo
        Std = sqrt(max(Std-E**2, 0d0))

    end function

end module
