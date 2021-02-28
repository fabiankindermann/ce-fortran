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

    ! number of subintervals
    integer, parameter :: n = 10

    ! old market prices
    real*8, parameter :: p_1 = 3d0

    ! new market price
    real*8, parameter :: p_2 = 1d0

    ! declaration of variables
    real*8 :: a, b, h, p
    real*8 :: d(0:n), w_trapez(0:n), f_trapez(0:n)
    real*8 :: d_gauss(0:n), w_gauss(0:n), f_gauss(0:n)
    real*8 :: trapez, simpson, gauss, CS(3, 2)

contains

    ! thefunction to store the inverse demand function
    function p_func(d_in)

        real*8, intent(in) :: d_in
        real*8 :: p_func

        ! the inverse demand function
        p_func = ((d_in + 1d0)/2d0)**(-2d0)

    end function

end module
