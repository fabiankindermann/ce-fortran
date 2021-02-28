!##############################################################################
! PROGRAM SummingUp
!
! ## Use do loops to sum up integer numbers successively
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
program SummingUp

    implicit none
    integer :: i, j, sum

    ! read the variable
    write(*,'(a)')'Type an integer variable that is > 0'
    read(*,*)i

    ! check whether i > 0
    if(i > 0)then

        ! initialize sum at 0
        sum = 0

        ! do the summing up
        do j = 1, i
            sum = sum + j
        enddo

        ! write the result to the console
        write(*,'(a,i3,a,i10)')'The sum of 1 to ',i,' is ',sum

    else

        ! write the error message
        write(*,'(a)')'Error: i should be greater than 0'
    endif

end program
