!##############################################################################
! PROGRAM Intertemporal
!
! ## Solve the intertemporal household optimization problem
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog02_02m.f90"

program Intertemporal

    use globals
    use toolbox

    implicit none

    ! declaration of variables
    real*8 :: x_root(3), x_minimize, a, b, fret
    logical :: check

    ! initialize x_root
    x_root = 0.1d0

    ! call subroutine fzero
    call fzero(x_root, foc, check)

    ! check whether fzero was executed correctly
    if(check)stop 'Error: fzero did not converge'

    ! print output
    write(*,'(a)') 'Result with fzero:'
    write(*,'(a)') '------------------'
    write(*,'(a,f10.6)') 'c_1    : ', x_root(1)
    write(*,'(a,f10.6)') 'c_2    : ', x_root(2)
    write(*,'(a,f10.6)') 'lambda : ', x_root(3)

    ! initialize interval for fminsearch and x_minimize
    a = 0d0
    b = w
    x_minimize = w/2d0

    ! call subroutine fminsearch
    call fminsearch(x_minimize, fret, a, b, utility)

    ! print output
    write(*,'(/a)') 'Result with fminsearch:'
    write(*,'(a)') '-----------------------'
    write(*,'(a,f10.6)') 'c_1    : ', x_minimize
    write(*,'(a,f10.6)') 'c_2    : ', (w - x_minimize)*(1d0+r)

end program
