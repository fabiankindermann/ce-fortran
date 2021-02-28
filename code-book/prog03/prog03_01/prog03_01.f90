!##############################################################################
! PROGRAM planner
!
! ## The social planner solution to the static general equilibrium model
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
include "prog03_01m.f90"

program planner

    use globals
    use toolbox

    implicit none
    real*8 :: x(2), L(2), K(2), Y(2), fret

    ! initial guess
    x(:) = 5d0

    ! minimization routine
    call fminsearch(x, fret, (/0d0, 0d0/), (/10d0, 20d0/), utility)

    ! solution
    K(1) = x(1)
    L(1) = x(2)
    K(2) = Kbar - K(1)
    L(2) = Lbar - L(1)
    Y    = L**beta*K**(1d0-beta)

    ! output
    write(*,'(/a)')'GOODS MARKET 1 :'
    write(*,'(4(a,f6.2,2x))')' X1 =',Y(1),' Y1 =',Y(1)

    write(*,'(/a)')'GOODS MARKET 2 :'
    write(*,'(4(a,f6.2,2x))')' X2 =',Y(2),' Y2 =',Y(2)

    write(*,'(/a)')'LABOR MARKET :'
    write(*,'(4(a,f6.2,2x))')' L1 =',L(1),' L2 =',L(2),' L  =',Lbar

    write(*,'(/a)')'CAPITAL MARKET :'
    write(*,'(4(a,f6.2,2x))')' K1 =',K(1),' K2 =',K(2),' K  =',Kbar

end program
