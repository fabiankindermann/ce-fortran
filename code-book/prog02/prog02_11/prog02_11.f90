!##############################################################################
! PROGRAM brentmin
!
! ## Finding the minimum of a function with simple constraints
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
include "prog02_11m.f90"

program brentmin

    use globals
    use toolbox

    implicit none
    real*8 :: x, f, a, b

    ! initial interval and function values
    a = 0d0
    b = (W-p(1)*0.01d0)/p(2)

    ! set starting point
    x = (a+b)/2d0

    ! call minimizing routine
    call fminsearch(x, f, a, b, utility)

    ! output
    write(*,'(/a,f12.7)')' x_1 = ',(W-p(2)*x)/p(1)
    write(*,'(a,f12.7)')' x_2 = ',x
    write(*,'(a,f12.7)')' u   = ',-f

end program
