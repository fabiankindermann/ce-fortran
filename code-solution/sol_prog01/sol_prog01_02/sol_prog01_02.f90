!##############################################################################
! PROGRAM RealAndStar8
!
! ## Compare type real with real*8
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
program RealAndStar8

    implicit none

    ! declaration of variables
    real :: sum1
    real*8 :: sum2

    ! calculation of the sums
    sum1 = 55555553 + 10000001
    sum2 = 55555553 + 10000001

    ! print output
    write(*,'(a, f15.2)')'Using real   =', sum1
    write(*,'(a, f15.2)')'Using real*8 =', sum2

end program
