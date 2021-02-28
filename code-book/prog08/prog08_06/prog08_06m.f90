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
    real*8, parameter :: beta = 0.95d0
    real*8, parameter :: a0 = 100d0

    ! numerical parameters
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! time path of consumption and resource
    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT), a_t(0:TT)

    ! policy function
    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NA), coeff_c(NA+3)
    real*8 :: con_lev, x_in
    logical :: check

    ! variables to communicate with function
    integer :: ia_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, cplus

        ! calculate right hand side of foc
        cplus = spline_eval(x_in, coeff_c, 0d0, a0)

        ! get foc
        foc = a(ia_com) - x_in - beta**(-gamma)*cplus

    end function

end module
