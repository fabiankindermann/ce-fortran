!##############################################################################
! PROGRAM TradingStrategies
!
! ## Options as part of trading strategies
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
program TradingStrategies

    use toolbox

    implicit none

    real*8, parameter :: Del_TT_3 = 21d0/250d0   ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 100d0             ! initial stock price
    real*8, parameter :: r = 0.1d0               ! annual interest rate
    real*8, parameter :: sigma = 0.4d0           ! standard deviation of stock returns
    real*8, parameter :: KK = 100d0              ! strike price straddle
    real*8, parameter :: KK_1_bfly = 95d0        ! strike price for butterfly spread
    real*8, parameter :: KK_1_strg = 105d0       ! strike price for strangle
    real*8, parameter :: KK_2_bfly = 105d0       ! strike price for butterfly spread
    real*8, parameter :: KK_2_strg = 95d0        ! strike price for strangle
    integer, parameter :: TT = 21                ! number of sub-periods (total)
    integer, parameter :: nplot = 20             ! number of gridpoints for plot

    real*8 :: c_E, p_E, c_E1, p_E1, c_E2, p_E2, KK_temp, d_1, d_2
    real*8 :: pi_strd(2), pi_strg(2), pi_bfly(4)
    real*8 :: p_strd, p_strg, p_bfly
    real*8 :: S_plot(0:nplot), pi_plot(0:nplot)
    integer :: in

    ! compute fair value of a straddle
    d_1 = (log(S_0/KK) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    c_E = S_0*normalCDF(d_1, 0d0, 1d0) - KK*exp(-r*Del_TT_3)*normalCDF(d_2, 0d0, 1d0)
    p_E = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK*exp(-r*Del_TT_3)*normalCDF(-d_2, 0d0, 1d0)
    p_strd = c_E + p_E

    ! plot output
    do in = 0, nplot

        ! set up x-axis
        S_plot(in) = 90d0 + dble(in)

        ! compute payoff
        pi_strd(1) = max(S_plot(in) - KK, 0d0)
        pi_strd(2) = max(KK - S_plot(in), 0d0)

        ! compute payoff for straddle
        pi_plot(in) = pi_strd(1) + pi_strd(2)

    enddo

    ! plot output
    call plot(S_plot, pi_plot, legend='European Straddle')
    call execplot()

    ! compute fair value of a strangle
    d_1 = (log(S_0/KK_1_strg) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    c_E = S_0*normalCDF(d_1, 0d0, 1d0) - KK_1_strg*exp(-r*Del_TT_3)*normalCDF(d_2, 0d0, 1d0)

    d_1 = (log(S_0/KK_2_strg) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    p_E = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK_2_strg*exp(-r*Del_TT_3)*normalCDF(-d_2, 0d0, 1d0)
    p_strg = c_E + p_E

    ! plot output
    do in = 0, nplot

        ! set up x-axis
        S_plot(in) = 90d0 + dble(in)

        ! compute payoff
        pi_strg(1) = max(S_plot(in) - KK_1_strg, 0d0)
        pi_strg(2) = max(KK_2_strg - S_plot(in), 0d0)

        ! compute payoff for strangle
        pi_plot(in) = pi_strg(1) + pi_strg(2)

    enddo

    ! plot output
    call plot(S_plot, pi_plot, legend='European Strangle')
    call execplot()

    ! compute fair value of a butterfly spread
    d_1 = (log(S_0/KK_1_bfly) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    c_E1 = S_0*normalCDF(d_1, 0d0, 1d0) - KK_1_bfly*exp(-r*Del_TT_3)*normalCDF(d_2, 0d0, 1d0)

    d_1 = (log(S_0/KK_2_bfly) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    p_E1 = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK_2_bfly*exp(-r*Del_TT_3)*normalCDF(-d_2, 0d0, 1d0)

    KK_temp = 0.5d0*(KK_1_bfly + KK_2_bfly)
    d_1 = (log(S_0/KK_temp) + (r+sigma**2/2d0)*Del_TT_3)/sigma/sqrt(Del_TT_3)
    d_2 = d_1 - sigma*sqrt(Del_TT_3)
    c_E2 = S_0*normalCDF(d_1, 0d0, 1d0) - KK_temp*exp(-r*Del_TT_3)*normalCDF(d_2, 0d0, 1d0)
    p_E2 = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK_temp*exp(-r*Del_TT_3)*normalCDF(-d_2, 0d0, 1d0)
    p_bfly = c_E1 + p_E1 - c_E2 - p_E2

    ! plot output
    do in = 0, nplot

        ! set up x-axis
        S_plot(in) = 90d0 + dble(in)

        ! compute payoff
        pi_bfly(1) = max(S_plot(in) - KK_1_bfly, 0d0)
        pi_bfly(2) = max(KK_2_bfly - S_plot(in), 0d0)
        pi_bfly(3) = max(S_plot(in) - KK_temp, 0d0)
        pi_bfly(4) = max(KK_temp - S_plot(in), 0d0)

        ! compute straddle payoff
        pi_plot(in) = pi_bfly(1) + pi_bfly(2) - pi_bfly(3) - pi_bfly(4)

    enddo

    ! plot output
    call plot(S_plot, pi_plot, legend='Butterfly Spread')
    call execplot(ylim=(/0d0, 12d0/))

    write(*,'(a)')'TRADING STRATEGY PRICING USING BLACK-SCHOLES FORMULA'
    write(*,'(a, f12.5)')'Straddle Price         = ', p_strd
    write(*,'(a, f12.5)')'Strangle Price         = ', p_strg
    write(*,'(a, f12.5)')'Butterfly Spread Price = ', p_bfly

end program
