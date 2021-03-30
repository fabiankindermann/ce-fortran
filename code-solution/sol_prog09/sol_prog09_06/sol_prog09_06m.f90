!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    implicit none

    ! number of transitory shock process values
    integer, parameter :: NS = 15

    ! risk process
    real*8, parameter :: sigma = 1d0
    real*8, parameter :: rho = 0.5d0

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 100

    ! counter variables
    integer :: iter, is, is_p

    ! the risk process
    real*8 :: pi(NS, NS), eta(NS)

    ! distribution variables
    real*8 :: phi(NS), phi_new(NS), phi_cdf(NS)
    real*8 :: con_lev

end module
