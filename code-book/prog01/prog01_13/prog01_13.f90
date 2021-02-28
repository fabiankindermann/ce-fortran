!##############################################################################
! PROGRAM Functions
!
! ## Store reusable code in functions with output argument
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
program Functions

    implicit none
    real*8 :: a, b, res

    a = 3d0
    b = 5d0

    ! call function
    res =  addIt(a,b)

    ! output
    write(*,'(2(f8.2,a),f8.2)')a,' + ',b,' = ',res

contains


    function addIt(a, b)

        implicit none

        ! input arguments
        real*8, intent(in) :: a, b

        ! function value
        real*8 :: addIt

        ! executable code
        addIt = a + b

    end function

end program
