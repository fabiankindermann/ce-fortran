!##############################################################################
! PROGRAM ImpliedVolatility
!
! ## The implied volatility of an european option for a given price cE
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog04_02m.f90"

program ImpliedVolatility

    use globals
    use toolbox

    implicit none
    real*8 :: sigma
    logical :: check

    ! set up initial guess
    sigma = 0.5d0

    ! call subroutine fzero
    call fzero(sigma, DIFF, check)

    ! check for convergence
    if(check)stop 'Error: fzero did not converge'

    ! print output to the screen
    write(*,'(a,f10.1)')'Implied Volatility (in %): ', sigma*100

end program
