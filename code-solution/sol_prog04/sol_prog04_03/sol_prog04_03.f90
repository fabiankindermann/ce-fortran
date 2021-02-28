!##############################################################################
! PROGRAM MCEuropean
!
! ## Option pricing with Monte Carlo simulations - European options
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
program MCEuropean

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    integer, parameter :: JJ = 20000             ! maximum number of simulation paths

    real*8 :: z, d_1, d_2
    real*8 :: S_T, pi_c(JJ), pi_p(JJ)
    real*8 :: c_E(JJ), p_E(JJ), c_EMC(JJ), p_EMC(JJ), npaths(JJ)
    integer :: ij

    ! use Black-Scholes formula for call and put price
    d_1 = (log(S_0/KK) + (r+sigma**2/2d0)*Del_TT)/(sigma*sqrt(Del_TT))
    d_2 = d_1 - sigma*sqrt(Del_TT)
    c_E = S_0*normalCDF(d_1, 0d0, 1d0) - KK*exp(-r*Del_TT)*normalCDF(d_2, 0d0, 1d0)
    p_E = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK*exp(-r*Del_TT)*normalCDF(-d_2, 0d0, 1d0)

    ! set up monte-carlo algorithm
    do ij = 1, JJ

        ! simulate stock-prices at maturity
        call simulate_normal(z, 0d0, 1d0)
        S_T = S_0*exp((r-sigma**2/2d0)*Del_TT + sigma*sqrt(Del_TT)*z)

        ! calculate payoff at maturity
        pi_c(ij) = max(S_T - KK, 0d0)
        pi_p(ij) = max(KK - S_T, 0d0)

    enddo

    ! use different number of paths to compute premium
    do ij = 1, JJ
        c_EMC(ij) = exp(-r*Del_TT)*sum(pi_c(1:ij))/dble(ij)
        p_EMC(ij) = exp(-r*Del_TT)*sum(pi_p(1:ij))/dble(ij)
        npaths(ij) = dble(ij)
    enddo

    ! plot option price as a function of number of paths used
    call plot(npaths(1:10000), c_EMC(1:10000), legend='European Call Option (MC)')
    call plot(npaths(1:10000), c_E(1:10000), legend='European Call Option')
    call plot(npaths(1:10000), p_EMC(1:10000), legend='European Put Option (MC)')
    call plot(npaths(1:10000), p_E(1:10000), legend='European Put Option')
    call execplot(xlabel='Number of simulated paths', ylabel='Simulated option price')

    write(*,'(a)')'OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Call Price = ', c_EMC(JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_EMC(JJ)
    write(*,*)
    write(*,'(a)')'OPTION PRICING USING BLACK-SCHOLES'
    write(*,'(a, f12.5)')'Call Price = ', c_E(JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_E(JJ)

end program
