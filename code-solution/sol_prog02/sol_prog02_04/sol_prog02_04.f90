!##############################################################################
! PROGRAM Intertemporal2dim
!
! ## Solve 3 period intertemporal optimization with 2 dimensional fminsearch
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
include "sol_prog02_04m.f90"

program Intertemporal2dim

    use globals
    use toolbox

    implicit none

    ! declaration of variables
    real*8 :: x_in(2), a(2), b(2), fret

    ! initialize interval for fminsearch and x_minimize
    a = (/0d0, 0d0/)
    b = (/(w + w*(1d0+r)), (w*(1d0+r) + w*(1d0+r)**2d0)/)
    x_in = w/2d0

    ! call subroutine fminsearch
    call fminsearch(x_in, fret, a, b, utility)

    ! print output
    write(*,'(/a)') 'Result with fminsearch:'
    write(*,'(a)') '-----------------------'
    write(*,'(a,f10.6)') 'c_1    : ', w + w/(1d0+r) - x_in(1)/(1d0+r) &
        - x_in(2)/(1d0+r)**2d0
    write(*,'(a,f10.6)') 'c_2    : ', x_in(1)
    write(*,'(a,f10.6)') 'c_3    : ', x_in(2)

end program
