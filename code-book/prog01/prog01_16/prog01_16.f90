!##############################################################################
! PROGRAM CalcUtil
!
! ## Calculate utility of household for different values of c1 and c2
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
include "prog01_16m.f90"

program CalcUtil

    use Globals
    use UtilFunc

    implicit none
    real*8 :: c1, c2, util
    integer :: j

    ! initialize parameters
    beta = 0.9d0
    eta = 2d0

    ! calculate utility for different consumption pairs
    ! between 0.3 and 0.7
    do j = 0, 20
        c1 = 0.3d0 + (0.7d0-0.3d0)/20*dble(j)
        c2 = 1d0-c1
        util = utility(c1, c2)
        write(*,'(3f10.4)')c1, c2, util
    enddo

end program
