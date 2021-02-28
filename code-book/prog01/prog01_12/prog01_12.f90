!##############################################################################
! PROGRAM Subroutines
!
! ## Store reusable code in subroutines
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
program Subroutines

    implicit none
    real*8 :: a, b, c, d

    a = 3d0
    b = 5d0

    ! call subroutine
    call addIt(a, b)

    ! redefine values
    c = 10d0
    d = 2d0

    ! call subroutine again
    call addIt(c, d)

! separates main program code from subroutine and functions
contains


    subroutine addIt(a, b)

        implicit none

        ! input arguments
        real*8, intent(in) :: a, b

        ! other variables
        real*8 :: c

        ! executable code
        c = a + b
        write(*,'(2(f8.2,a),f8.2)')a,' + ',b,' = ',c

    end subroutine

end program
