!##############################################################################
! PROGRAM IfStatements
!
! ## How to use if statements to control program flow
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
program IfStatements

    implicit none
    integer :: a

    ! initialize a
    a = 1

    ! check for the size of a
    if(a < 1)then
        write(*,'(a)')'condition 1 is true'
    elseif(a < 2)then
        write(*,'(a)')'condition 2 is true'
    elseif(a < 3)then
        write(*,'(a)')'condition 3 is true'
    else
        write(*,'(a)')'no condition is true'
    endif

end program
