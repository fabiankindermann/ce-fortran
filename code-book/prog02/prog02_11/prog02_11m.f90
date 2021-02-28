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
    real*8, parameter :: p(2) = (/1d0, 2d0/)
    real*8, parameter :: W = 1d0

contains


    ! the utility function
    function utility(x)

        implicit none
        real*8, intent(in) :: x
        real*8 :: utility

        utility = -(((W-p(2)*x)/p(1))**0.4d0+(1d0+x)**0.5d0)

    end function

end module globals
