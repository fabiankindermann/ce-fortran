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

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: KK = 25d0               ! strike price
    real*8, parameter :: c_E = 0.63d0            ! option value

contains

    function DIFF(sigma)

        implicit none
        real*8, intent(in) :: sigma
        real*8 :: DIFF, d_1, d_2

        d_1 = (log(S_0/KK) + (r+sigma**2/2d0)*Del_TT)/sigma/sqrt(Del_TT)
        d_2 = d_1 - sigma*sqrt(Del_TT)
        DIFF = S_0*normalCDF(d_1, 0d0, 1d0) - &
               KK*exp(-r*Del_TT)*normalCDF(d_2, 0d0, 1d0) - c_E

    end function

end module
