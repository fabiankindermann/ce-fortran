!##############################################################################
! PROGRAM BinaryOption
!
! ## Use closed-form solution to price a Binary option
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
program BinaryOption

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    real*8, parameter :: q = 30                  ! payoff at maturity

    real*8 :: c_BI, p_BI, d_2

    ! compute price of Binary call and put
    d_2 = (log(S_0/KK) + (r-sigma**2/2d0)*Del_TT)/sigma/sqrt(Del_TT)
    c_BI = q*normalCDF(d_2, 0d0, 1d0)
    p_BI = q*normalCDF(-d_2, 0d0, 1d0)

    ! write output
    write(*,'(a,f10.2)')'Binary Call: ', c_BI
    write(*,'(a,f10.2)')'Binary Put:  ', p_BI

end program
