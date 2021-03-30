!##############################################################################
! PROGRAM SingleAndDoublePrecision
!
! ## Compare single and double precision with each other
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program SingleAndDoublePrecision

    implicit none

    ! declaration of variables
    real*8 :: rvar1(3), rvar2(3), x(2), y(2)
    integer :: i

    ! assign 10**9 to rvar1
    rvar1(1) = 10**9
    rvar1(2) = 10d0**9
    rvar1(3) = 10**9d0

    ! assign 10**10 to rvar2
    !rvar2(1) = 10**10
    rvar2(2) = 10d0**10
    rvar2(3) = 10**10d0

    ! print output
    write(*,'(a)')'                  Without d0          With base d0          With expo d0'
    write(*,'(a, 3f22.2)')'Exp  9', (rvar1(i), i = 1, 3)
    write(*,'(a, 3f22.2)')'Exp 10', (rvar2(i), i = 1, 3)

    ! assign 0.000000000003 to x with single and double precision
    x(1) = 0.000000000003d0
    x(2) = 0.000000000003

    write(*,'(//a, 2f30.25)')'Precision   ', (x(i), i = 1, 2)

    ! assign 3.1415926535 to y with single and double precision
    y(1) = 3.1415926535d0
    y(2) = 3.1415926535

    write(*,'(//a, 2f15.12)')'Precision     ', (y(i), i = 1, 2)

end program
