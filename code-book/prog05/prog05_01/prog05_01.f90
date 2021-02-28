!##############################################################################
! PROGRAM household1
!
! ## The life cycle model without uncertainty
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
include "prog05_01m.f90"

program household1

    use globals
    use toolbox

    implicit none
    real*8 :: x(2)
    real*8 :: fret, low(2), up(2)

    ! lower and upper border and initial guess
    low = (/0d0, 0d0/)
    up = (/w, R*w+w/)
    x = up/2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC    SAV'
    write(*,'(i4,4f7.2/)')1,c(1),w,w,a(2)
    write(*,'(i4,4f7.2/)')2,c(2),w,w+R*a(2),a(3)
    write(*,'(i4,4f7.2)')3,c(3),0d0,R*a(3),0d0

end program
