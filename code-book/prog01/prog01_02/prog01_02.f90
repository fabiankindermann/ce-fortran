!##############################################################################
! PROGRAM Hello
!
! ## Read and write numbers from/to the console
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
program Hello

    ! declaration of variables
    implicit none
    character(len=50) :: input

    ! executable code
    write(*,*)'Please type your name:'
    read(*,*)input
    write(*,*)'Hello ',input

end program
