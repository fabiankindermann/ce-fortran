!##############################################################################
! MODULE globals
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
module globals

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: beta = 0.95d0
    real*8, parameter :: a0 = 100d0

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ia, ia_p, iter

    ! time path of consumption and resource
    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT), a_t(0:TT)
    integer :: ia_t(0:TT)

    ! value and policy function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA), V(0:NA)

    ! variables to numerically determine value and policy function
    integer :: ia_opt(0:NA)
    real*8 :: V_new(0:NA), u_temp(0:NA)
    real*8 :: con_lev, cons

end module
