!##############################################################################
! PROGRAM VarDec
!
! ## Shows how to declare variables
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
program VarDec

    ! declaration of variables
    implicit none
    logical :: logic
    integer :: a, b
    real*8 :: x, y1
    character :: one_char
    character(len=20) :: long_char

    real*8, parameter :: pi = 3.14d0
    integer, parameter :: n = 56

end program
