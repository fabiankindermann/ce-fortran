!##############################################################################
! PROGRAM Logicals
!
! ## Evaluate logical expression
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
program Logicals

    implicit none

    ! declare variables
    integer :: x, y, z

    ! assign values
    x = 4
    y = 6
    z = 8

    ! print output
    write(*,*)(x >= 3 .and. y <= 4 .and. z == 5 .or. x <= y .and. y < z)

end program
