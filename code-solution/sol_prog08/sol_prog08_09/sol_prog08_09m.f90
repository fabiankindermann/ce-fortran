!##############################################################################
! MODULE globals
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
module globals

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: beta = 0.975d0
    real*8, parameter :: r = 0.02d0
    real*8, parameter :: w = 1d0
    real*8, parameter :: a0 = 0d0
    real*8, parameter :: a_borrow = 5d0

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ia, iter

    ! time path of consumption and resource
    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT), a_t(0:TT)

    ! policy function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA)
    real*8 :: a_l, a_u

    ! variables to numerically determine policy function
    integer :: ia_low
    real*8 :: a_endog(-1:NA), c_endog(-1:NA)
    real*8 :: c_new(0:NA), coeff_c(NA+3)
    real*8 :: con_lev

end module
