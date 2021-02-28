!##############################################################################
! PROGRAM Arrays
!
! ## Use do-loops to fill up arrays with values
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
program Arrays

    implicit none
    real*8 :: x(0:10), y(0:10)
    integer :: j

    ! initialize x and calculate y
    do j = 0, 10
        x(j) = 1d0/10d0*dble(j)
        y(j) = exp(x(j))
    enddo

    ! output table of values
    write(*,'(a)')'         X         Y'
    do j = 0, 10
        write(*,'(2f10.3)')x(j), y(j)
    enddo

end program
