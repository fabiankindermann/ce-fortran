!##############################################################################
! PROGRAM Chooser
!
! ## Analytical option pricing and Monte Carlo simulations - Chooser options
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program Chooser

    use toolbox

    implicit none
    real*8, parameter :: Del_TT_1 = 31d0/250d0   ! exercise date (in annualized values)
    real*8, parameter :: Del_TT_2 = 62d0/250d0   ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    integer, parameter :: TT_1 = 31              ! number of sub-periods
    integer, parameter :: TT_2 = 62              ! number of sub-periods
    integer, parameter :: JJ = 20000             ! maximum number of simulation paths

    real*8 :: Dis_KK, S(0:TT_1), z(TT_1), del_t, pi(JJ)
    real*8 :: p_CH(JJ), p_CHMC(JJ), npaths(JJ), c_E, p_E, d_1, d_2
    integer :: ij, it

    ! compute parameters for Black-Scholes formula with maturity T2
    d_1 = (log(S_0/KK) + (r+sigma**2/2d0)*Del_TT_2)/sigma/sqrt(Del_TT_2)
    d_2 = d_1 - sigma*sqrt(Del_TT_2)
    c_E = S_0*normalCDF(d_1, 0d0, 1d0) &
          - KK*exp(-r*Del_TT_2)*normalCDF(d_2, 0d0, 1d0)

    ! compute parameters for Black-Scholes formula with maturity T1
    Dis_KK = KK*exp(-r*(Del_TT_2-Del_TT_1))
    d_1 = (log(S_0/Dis_KK) + (r+sigma**2/2d0)*Del_TT_1)/sigma/sqrt(Del_TT_1)
    d_2 = d_1 - sigma*sqrt(Del_TT_1)
    p_E = -S_0*normalCDF(-d_1, 0d0, 1d0) &
          + KK*exp(-r*Del_TT_1)*normalCDF(-d_2, 0d0, 1d0)

    ! compute the price of a Chooser option
    p_CH = c_E + p_E

    ! Monte-Carlo simulation
    del_t = Del_TT_1/dble(TT_1)

    do ij = 1, JJ

        ! generate a price path
        call simulate_normal(z, 0d0, 1d0)
        S(0) = S_0
        ! simulate paths of stock prices until T1
        do it = 1, TT_1
            S(it) = S(it - 1)*exp((r-sigma**2/2d0)*del_t &
                                  + sigma*sqrt(del_t)*z(it))
        enddo
        d_1 = (log(S(TT_1)/KK) + (r+sigma**2/2d0)*(Del_TT_2-Del_TT_1)) &
              /sigma/sqrt(Del_TT_2-Del_TT_1)
        d_2 = d_1 - sigma*sqrt(Del_TT_2-Del_TT_1)
        c_E = S(TT_1)*normalCDF(d_1, 0d0, 1d0) &
              - KK*exp(-r*(Del_TT_2-Del_TT_1))*normalCDF(d_2, 0d0, 1d0)
        p_E = -S(TT_1)*normalCDF(-d_1, 0d0, 1d0) &
              + KK*exp(-r*(Del_TT_2-Del_TT_1))*normalCDF(-d_2, 0d0, 1d0)

        ! compute max price at T1
        pi(ij) = max(c_E, p_E)

    enddo

    ! compute price for lookback options from simulated paths
    do ij = 1, JJ
        p_CHMC(ij) = exp(-r*Del_TT_1)*sum(pi(1:ij))/dble(ij)
        npaths(ij) = dble(ij)
    enddo

    ! plot Chooser option price as a function of number of paths used
    call plot(npaths(1:10000), p_CHMC(1:10000), legend='Chooser Option (MC)')
    call plot(npaths(1:10000), p_CH(1:10000), legend='Chooser Option')
    call execplot()

    write(*,*)p_CHMC(10000)

    write(*,'(a)')'CHOOSER OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Price = ', p_CHMC(JJ)
    write(*,*)
    write(*,'(a)')'OPTION PRICING USING ANALYCTICAL SOLUTION'
    write(*,'(a, f12.5)')'Price = ', p_CH(JJ)

end program
