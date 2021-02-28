!##############################################################################
! MODULE Volume
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
module Volume

    implicit none

    ! module parameter
    real*8, parameter :: pi = 3.14159265358d0

! separates variable declarations from subroutine and functions
contains


    ! for calculating volume of a sphere
    function vol(r)

        implicit none

        ! input and output variables
        real*8, intent(in) :: r
        real*8 :: vol

        ! calculation
        vol = 4d0/3d0*r**3*pi

    end function

end module
