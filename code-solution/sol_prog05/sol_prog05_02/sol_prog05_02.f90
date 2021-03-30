!##############################################################################
! PROGRAM housing
!
! ## The life cycle model with uncertainty on health expenditures
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog05_02m.f90"

program healthexpend

    use globals
    use toolbox

    implicit none
    real*8 :: x(2)
    real*8 :: fret, low(2), up(2)

    ! discretize h
    if(sigma /= 0d0)then
        call grid_Cons_Equi(hc, mu-sigma, mu+sigma)
    else
        hc = mu
    endif
    weight_h = 1d0/NH

    ! lower and upper border and initial guess
    low = (/0d0, 0d0/)
    up = (/w, R*w + w/)
    x = up/2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC    SAV   EH'
    write(*,'(i4,4f7.2/)')1,c(1),w,w,a(2)
    write(*,'(i4,4f7.2/)')2,c(2),w,w+R*a(2),a(3)
    write(*,'(i4,5f7.2,a)')3,E(c3(:)),0d0,R*a(3),0d0,E(hc(:)),' (MEAN)'
    write(*,'(4x,5f7.2,a/)')Std(c3(:)),0d0,0d0,0d0,Std(hc(:)),' (STD)'

end program
