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
    real*8, parameter :: beta = 0.99d0
    real*8, parameter :: alpha_c = 0.5d0
    real*8, parameter :: alpha_i = 0.5d0
    real*8, parameter :: delta = 0.075d0
    real*8, parameter :: k0 = 10d0

    ! numerical parameters
    real*8, parameter :: k_l = 5d0
    real*8, parameter :: k_u = 100d0
    real*8, parameter :: sig = 1d-6
    integer, parameter :: itermax = 2000

    ! counter variables
    integer :: it, ik, iter

    ! time path of consumption and capital
    integer, parameter :: TT = 360
    real*8 :: c_t(0:TT), i_t(0:TT), k_t(0:TT), kc_t(0:TT), ki_t(0:TT)

    ! policy function
    integer, parameter :: NK = 100
    real*8 :: k(0:NK), kc(0:NK)

    ! variables to numerically determine policy function
    real*8 :: kc_new(0:NK), coeff_kc(NK+3)
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
        real*8 :: foc, kplus, kcplus, kc

        ! limit current consumption capital
        kc = max(x_in, 1d-10)
        kc = min(kc, k_com - 1d-10)

        ! future capital
        kplus = (1d0-delta)*k_com + (k_com-kc)**alpha_i

        ! calculate future consumption
        kcplus = spline_eval(kplus, coeff_kc, k_l, k_u)
        kcplus = max(kcplus, 1d-10)
        kcplus = min(kcplus, kplus - 1d-10)

        ! get first order condition
        foc = kc - (beta*(1d0-delta+alpha_i*(kplus - kcplus) &
                          **(alpha_i-1d0))*(kc/kcplus)**(1d0-alpha_c) &
                       *((kplus - kcplus)/(k_com - kc))**(1d0-alpha_i)) &
                       **(-gamma/alpha_c)*kcplus

    end function

end module
