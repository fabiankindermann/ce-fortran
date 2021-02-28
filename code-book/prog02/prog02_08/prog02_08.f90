!##############################################################################
! PROGRAM oligopoly
!
! ## Multi-dimensional rootfinding with the toolbox routine fzero
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
include "prog02_08m.f90"

program oligopoly

    use globals
    use toolbox

    implicit none
    real*8 :: q(2)
    logical :: check

    ! initialize q
    q = 0.1d0

    ! find root
    call fzero(q, cournot, check)

    if(check)stop 'Error: fzero did not converge'

    ! output
    write(*,'(/a)')'            Output'
    write(*,'(a,f10.4)')'Firm 1: ',q(1)
    write(*,'(a,f10.4)')'Firm 2: ',q(2)
    write(*,'(/a,f10.4)')'Price : ',(q(1)+q(2))**(-1d0/1.6d0)

end program
