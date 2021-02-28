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
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha = 0.40d0
    real*8, parameter :: delta = 0.019d0

    ! government parameters (change accordingly to simulate Table 9.1)
    real*8, parameter :: tau_r = 0.00d0
    real*8, parameter :: gy = 0.00d0
    real*8, parameter :: by = 0.00d0

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 100d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, iter

    ! government consumption
    real*8 :: gbar

    ! time path of consumption and capital
    integer, parameter :: TT = 3600
    real*8 :: c_t(0:TT), i_t(0:TT), k_t(0:TT), y_t(0:TT), w_t(0:TT), r_t(0:TT)
    real*8 :: tauw_t(0:TT), b_t(0:TT)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), c(0:NK)

    ! variables to numerically determine policy function
    real*8 :: c_new(0:NK), coeff_c(NK+3)
    real*8 :: con_lev, x_in
    logical :: check

    ! variables to communicate with function
    real*8 :: k_com

contains


    ! the first order condition
    function foc(x_in)

        use toolbox

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc, kplus, bplus, cplus

        ! future capital
        kplus = (1d0-delta)*k_com + k_com**alpha - x_in - gbar

        ! calculate future consumption
        cplus = spline_eval(kplus, coeff_c, k_l, k_u)

        ! get first order condition
        foc = x_in - (beta*(1d0+(1d0-tau_r)*(alpha*kplus**(alpha-1d0)-delta)))**(-gamma)*cplus

    end function

end module
