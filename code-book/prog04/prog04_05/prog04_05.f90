!##############################################################################
! PROGRAM MC_Option_Pricing
!
! ## Option pricing with Monte Carlo simulations - Asian options
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
program MC_Option_Pricing

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0                ! strike price
    integer, parameter :: TT = 200               ! number of sub-periods
    integer, parameter :: JJ = 20000             ! maximum number of simulation paths

    real*8 :: S(0:TT), z(TT), del_t
    real*8 :: S_bar, pi_c(JJ), pi_p(JJ)
    real*8 :: c_AS(JJ), p_AS(JJ), npaths(JJ)
    integer :: ij, it

    ! get length of sub-periods
    del_t = Del_TT/dble(TT)

    ! simulate price paths and corresponding profits
    do ij = 1, JJ

        ! generate a price path
        call simulate_normal(z, 0d0, 1d0)
        S(0) = S_0
        do it = 1, TT
            S(it) = S(it-1)*exp((r-sigma**2/2d0)*del_t + sigma*sqrt(del_t)*z(it))
        enddo

        ! get average stock price
        S_bar = 0
        do it = 1, TT
            S_bar = S_bar + S(it)
        enddo
        S_bar = S_bar/dble(TT)

        ! calculate associated payoff
        pi_c(ij) = max(S_bar - KK, 0d0)
        pi_p(ij) = max(KK - S_bar, 0d0)
    enddo

    ! use different number of paths to compute premium
    do ij = 1, JJ
        c_AS(ij) = exp(-r*Del_TT)*sum(pi_c(1:ij))/dble(ij)
        p_AS(ij) = exp(-r*Del_TT)*sum(pi_p(1:ij))/dble(ij)
        npaths(ij) = dble(ij)
    enddo

    ! plot option price as a function of number of paths used
    call plot(npaths(1:10000), c_AS(1:10000), legend='Asian Call Option')
    call plot(npaths(1:10000), p_AS(1:10000), legend='Asian Put Option')
    call execplot()

    write(*,'(a)')'OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Call Price = ', c_AS(JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_AS(JJ)

end program
