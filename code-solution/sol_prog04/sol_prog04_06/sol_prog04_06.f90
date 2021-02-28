!##############################################################################
! PROGRAM MCLookback
!
! ## Option pricing with Monte Carlo simulations - Lookback options
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
program MCLookback

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    integer, parameter :: TT = 200               ! number of sub-periods
    integer, parameter :: JJ = 20000             ! maximum number of simulation paths
    integer, parameter :: NN = 2                 ! number of different option types

    real*8 :: S(0:TT), z(TT), del_t
    real*8 :: pi_c(NN, JJ), pi_p(NN, JJ)
    real*8 :: c_LB(NN, JJ), p_LB(NN, JJ), npaths(JJ)
    integer :: ij, in, it

    ! get length of sub-periods
    del_t = Del_TT/dble(TT)

    do ij = 1, JJ

        ! generate a price path
        call simulate_normal(z, 0d0, 1d0)
        S(0) = S_0
        do it = 1, TT
            S(it) = S(it-1)*exp((r-sigma**2/2d0)*del_t + sigma*sqrt(del_t)*z(it))
        enddo

        ! calculate payoff at maturity for a floating strike option
        pi_c(1, ij) = max(S(TT) - minval(S(:)), 0d0)
        pi_p(1, ij) = max(maxval(S(:)) - S(TT), 0d0)

        ! calculate payoff at maturity for a fixed strike option
        pi_c(2, ij) = max(maxval(S(:)) - KK, 0d0)
        pi_p(2, ij) = max(KK - minval(S(:)), 0d0)
    enddo

    ! compute price for lookback options from simulated paths
    do ij = 1, JJ
        do in = 1, NN
            c_LB(in, ij) = exp(-r*Del_TT)*sum(pi_c(in, 1:ij))/dble(ij)
            p_LB(in, ij) = exp(-r*Del_TT)*sum(pi_p(in, 1:ij))/dble(ij)
            npaths(ij) = dble(ij)
        enddo
    enddo

    ! plot floating strike option price as a function of number of paths used
    call plot(npaths(1:10000), c_LB(1, 1:10000), legend='Floating Strike Call Option')
    call plot(npaths(1:10000), p_LB(1, 1:10000), legend='Floating Strike Put Option')
    call execplot()

    ! plot fixed strike option price as a function of number of paths used
    call plot(npaths(1:10000), c_LB(2, 1:10000), legend='Fixed Strike Call Option')
    call plot(npaths(1:10000), p_LB(2, 1:10000), legend='Fixed Strike Put Option')
    call execplot()

    write(*,'(a)')'FLOATING STRIKE OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Call Price = ', c_LB(1, JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_LB(1, JJ)
    write(*,*)
    write(*,'(a)')'FIXED STRIKE OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Call Price = ', c_LB(2, JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_LB(2, JJ)

end program
