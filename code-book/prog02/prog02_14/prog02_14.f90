!##############################################################################
! PROGRAM agriculture
!
! ## Solve the agricultural problem with discretized normal distribution
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
include "prog02_14m.f90"

program agriculture

    use globals
    use toolbox

    implicit none
    real*8 :: A, Ep, Varp
    logical :: check

    ! discretize y
    call normal_discrete(y, w, mu, sig2)

    ! initialize variables
    A = 1d0

    ! get optimum
    call fzero(A, market, check)

    ! get expectation and variance of price
    Ep = sum(w*max(3d0-2d0*A*y, minp))
    Varp = sum(w*(max(3d0-2d0*A*y, minp) - Ep)**2)

    ! write output
    write(*,'(a,f10.3)')' A      = ',A
    write(*,'(a,f10.3)')' E(p)   = ',Ep
    write(*,'(a,f10.3)')' Var(p) = ',Varp

end program
