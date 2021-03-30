!##############################################################################
! PROGRAM housing
!
! ## The life cycle model with a durable good
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog05_01m.f90"

program housing

    use globals
    use toolbox

    implicit none
    real*8 :: x(3)
    real*8 :: fret, low(3), up(3)

    ! lower and upper border and initial guess
    low = (/0d0, -a_low, -a_low/)
    up = (/w + a_low, w, R*w + w/)
    x = up/3d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! output
    write(*,'(/a/)')' AGE   CONS   DCONS  WAGE    INC    SAV   UTIL'
    write(*,'(i4,5f7.2/)')1,c(1),ah,w,w,a(2)
    write(*,'(i4,5f7.2/)')2,c(2),ah,w,w+R*a(2)-delh*ah,a(3)
    write(*,'(i4,6f7.2)')3,c(3),ah,0d0,R*a(3)+(1d0-delh+r)*ah,0d0,-fret

end program
