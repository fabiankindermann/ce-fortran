!##############################################################################
! PROGRAM MCBarrier
!
! ## Option pricing with Monte Carlo simulations - Barrier options
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program MCBarrier

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    real*8, parameter :: H = 24d0                ! barrier
    integer, parameter :: TT = 62                ! number of sub-periods
    integer, parameter :: TT_B = 10              ! number of sub-periods in which barrier is active
    integer, parameter :: JJ = 20000             ! maximum number of simulation paths
    integer, parameter :: NN = 3                 ! number of different option types

    real*8 :: S(0:TT), z(TT), del_t
    real*8 :: pi_c(NN, JJ), pi_p(NN, JJ)
    real*8 :: c_EB(NN, JJ), p_EB(NN, JJ), npaths(JJ)
    integer :: ij, in, it
    logical :: active(NN)

    ! get length of sub-periods
    del_t = Del_TT/dble(TT)

    ! simulate price paths and corresponding profits
    do ij = 1, JJ

        ! reset barrier
        active = .false.

        ! generate a price path
        call simulate_normal(z, 0d0, 1d0)
        S(0) = S_0
        do it = 1, TT
            S(it) = S(it-1)*exp((r-sigma**2/2d0)*del_t &
                                  + sigma*sqrt(del_t)*z(it))

            ! check wether barrier is hit at all
            if(S(it) < H)active(1) = .true.
            if(S(it) < H .and. it < TT_B)active(2) = .true.
            if(S(it) < H .and. it > TT - TT_B)active(3) = .true.

        enddo

        pi_C(:, ij) = 0d0
        pi_P(:, ij) = 0d0

        ! compute payoff at maturity conditional on whether barrier is hit
        if(active(1))then
            pi_C(1, ij) = max(S(TT) - KK, 0d0)
            pi_P(1, ij) = max(KK - S(TT), 0d0)
        endif

        if(active(2))then
            pi_C(2, ij) = max(S(TT) - KK, 0d0)
            pi_P(2, ij) = max(KK - S(TT), 0d0)
        endif

        if(active(3))then
            pi_C(3, ij) = max(S(TT) - KK, 0d0)
            pi_P(3, ij) = max(KK - S(TT), 0d0)
        endif
    enddo

    ! compute price for lookback options from simulated paths
    do ij = 1, JJ
        do in = 1, NN
            c_EB(in, ij) = exp(-r*Del_TT)*sum(pi_c(in, 1:ij))/dble(ij)
            p_EB(in, ij) = exp(-r*Del_TT)*sum(pi_p(in, 1:ij))/dble(ij)
            npaths(ij) = dble(ij)
        enddo
    enddo

    ! plot floating strike option price as a function of number of paths used
    call plot(npaths(1:10000), c_EB(1, 1:10000), legend='Barrier Call Option')
    call plot(npaths(1:10000), p_EB(1, 1:10000), legend='Barrier Put Option')
    call execplot(xlabel='Number of simulated paths', ylabel='Simulated option price')

    ! plot fixed strike option price as a function of number of paths used
    call plot(npaths(1:10000), c_EB(2, 1:10000), legend='Barrier Call Option (TT_B)')
    call plot(npaths(1:10000), p_EB(2, 1:10000), legend='Barrier Put Option (TT_B)')
    call execplot(xlabel='Number of simulated paths', ylabel='Simulated option price')

    ! plot fixed strike option price as a function of number of paths used
    call plot(npaths(1:10000), c_EB(3, 1:10000), legend='Barrier Call Option (TT - TT_B)')
    call plot(npaths(1:10000), p_EB(3, 1:10000), legend='Barrier Put Option (TT - TT_B)')
    call execplot(xlabel='Number of simulated paths', ylabel='Simulated option price')

    write(*,'(a)')'BARRIER OPTION PRICING USING MC SIMULATION'
    write(*,'(a, f12.5)')'Call Price = ', c_EB(1, JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_EB(1, JJ)
    write(*,*)
    write(*,'(a)')'BARRIER OPTION PRICING USING MC SIMULATION (TT_B)'
    write(*,'(a, f12.5)')'Call Price = ', c_EB(2, JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_EB(2, JJ)
    write(*,*)
    write(*,'(a)') 'BARRIER OPTION PRICING USING MC SIMULATION (TT - TT_B)'
    write(*,'(a, f12.5)')'Call Price = ', c_EB(3, JJ)
    write(*,'(a, f12.5)')'Put  Price = ', p_EB(3, JJ)

end program
