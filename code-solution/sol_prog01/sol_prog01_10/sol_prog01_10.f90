!##############################################################################
! PROGRAM SubFunc
!
! ## Calculate utilty for different values of consumption within subroutine
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
! #VC# VERSION: 1.0  (16 December 2019)
!
!##############################################################################
include "sol_prog01_10m.f90"

program SubFunc

    use globals

    implicit none

    ! declaration of variables
    integer :: i
    real*8 :: u(n), a, b

    ! initialize values for a and b
    a = 0.5d0
    b = 5d0

    ! call subroutine to compute utilities
    call utility_int(a, b, u)

    ! print output
    write(*,'(f12.5)') (u(i), i = 1, n)

end program SubFunc
