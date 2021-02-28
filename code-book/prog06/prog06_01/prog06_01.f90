!##############################################################################
! PROGRAM LR_OLG
!
! ## Long-run equilibria in the overlapping generations model
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
include "prog06_01m.f90"

program LR_OLG

    use globals
    use toolbox

    implicit none
    logical :: check
    integer :: j
    real*8 :: x(2)

    ! initialize labor supply, pension payments and tax rates
    LL = (2d0+n_p)/(1d0+n_p)
    taup = kappa/((2d0+n_p)*(1d0+n_p))
    tauc = 0d0
    tauw = 0d0
    taur = 0d0

    ! get initial guess
    x(:) = 0.7d0

    ! solve the steady state equation system
    call fzero(x, eqns, check)

    ! check whether the solution is valid
    if(check)stop 'No equilibirium found !'

    ! calculate household utility
    util = 0d0
    do j = 1, 3
        util = util + beta**(j-1)*c(j)**egam/egam
    enddo

    ! Output
    write(*,'(a/)')'  Steady state equilibrium'
    write(*,'(a)')'     c1     c2     c3      Y      w      r      U '
    write(*,'(7f7.2/)')c(:), YY, w, r, util

    write(*,'(a)')'     a2     a3      K '
    write(*,'(3f7.2/)')a(2), a(3), KK

    write(*,'(a)')'   tauw   taur   tauc   taup    pen      B'
    write(*,'(6f7.2/)')tauw, taur, tauc, taup, pen, BB

    write(*,'(a)')'      Y      C      G      I      DIFF'
    write(*,'(4f7.2,1f10.4)')YY,CC,GG,II,YY-CC-GG-II

end program
