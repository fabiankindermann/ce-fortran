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
    real*8, parameter :: mu = 1d0, sig2 = 0.1d0, minp = 1d0
    integer, parameter :: n = 10
    real*8 :: y(0:n), w(0:n)

contains


    function market(A)

        implicit none
        real*8, intent(in) :: A
        real*8 :: market
        real*8 :: Ep

        ! calculate expected price
        Ep = sum(w*max(3d0-2d0*A*y, minp))

        ! get equilibrium equation
        market = A - (0.5d0+0.5d0*Ep)

    end function

end module globals
