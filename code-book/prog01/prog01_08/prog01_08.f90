!##############################################################################
! PROGRAM DoLoops
!
! ## How to use do-loops to control program flows
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
program DoLoops

    implicit none
    integer :: j

    ! perform a do-loop for j = 1 to 10
    do j = 1, 10
        write(*,'(i3)')j
    enddo

    ! write a blank line
    write(*,*)

    ! perform a do-loop for j = 10 to 1
    do j = 10, 1, -1
        write(*,'(i3)')j
    enddo

    ! write a blank line
    write(*,*)

    ! alternative do-loop
    j = 1
    do
        write(*,'(i3)')j
        j = j + 1

        ! exit the do-loop
        if(j > 10)exit
    enddo

end program
