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
    integer, parameter :: TT = 5        ! number of time periods
    integer, parameter :: NN = 3        ! number of stocks
    real*8, parameter :: gamma = 10d0   ! risk aversion
    real*8, parameter :: r_f = 0.04     ! the risk free rate

    real*8 :: mu(NN), sig(NN, NN)
    real*8 :: x_l(NN), x_u(NN)

contains


    ! the utility function (risk-free)
    function utility_rf(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility_rf
        real*8 :: w(NN)

        w = (/x(1), x(2), x(3)/)

        utility_rf = -(dot_product(w, mu-r_f) - 0.5d0*gamma*dot_product(w, matmul(sig, w)))

    end function

    ! the utility function
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility
        real*8 :: w(NN)

        w = (/x(1), x(2), 1d0-x(1)-x(2)/)

        utility = -(dot_product(w, mu) - 0.5d0*gamma*dot_product(w, matmul(sig, w)))

    end function

end module
