!##############################################################################
! PROGRAM ArrayFunc
!
! ## Hand over arrays to functions (works the same with subroutines)
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
program ArrayFunc

    implicit none
    real*8 :: a(2), b(5), res(2,5)

    a(:) = (/3d0, 4d0/)
    b(:) = (/1d0, 2d0, 3d0, 4d0, 5d0/)

    ! call function
    res =  addIt(a,b)

    ! output
    write(*,'(5f8.2/5f8.2)')res(1,:),res(2,:)

contains


    function addIt(a, b)

        implicit none

        ! input arguments
        real*8, intent(in) :: a(2), b(:)

        ! function value
        real*8 :: addIt(2, size(b))

        ! local variables
        integer :: j, k

        ! executable code
        do j = 1, 2
            do k = 1, size(b)
                addIt(j,k) = a(j) + b(k)
            enddo
        enddo

    end function

end program
