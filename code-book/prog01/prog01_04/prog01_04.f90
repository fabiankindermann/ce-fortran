!##############################################################################
! PROGRAM ReadWrite
!
! ## Read and write number of different types from console
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
program ReadWrite

    ! declaration of variables
    implicit none
    integer :: a
    real*8 :: x

    ! executable code
    write(*,*)'Type an integer number:'
    read(*,*)a

    write(*,*)'Type a real number:'
    read(*,*)x

    write(*,*)a, x

end program
