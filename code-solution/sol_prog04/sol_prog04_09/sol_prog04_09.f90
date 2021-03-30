!##############################################################################
! PROGRAM GreeksEuropean
!
! ## Options as part of trading strategies
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog04_09m.f90"

program GreeksEuropean

    use globals

    implicit none
    integer :: in
    real*8 :: delta_c_E, delta_p_E, gamma_c_E, gamma_p_E, c_E, p_E
    real*8 :: vega_c_E, vega_p_E, theta_c_E, theta_p_E, rho_c_E, rho_p_E
    real*8 :: S_plot(0:nplot), p_plot(2, 0:nplot), delta_plot(2, 0:nplot)
    real*8 :: gamma_plot(2, 0:nplot), vega_plot(2, 0:nplot)
    real*8 :: theta_plot(2, 0:nplot), rho_plot(2, 0:nplot)

    ! compute fair value of the option and sensitivities
    call black_scholes(S_0, r, sigma, Del_TT, KK, c_E, p_E)
    call approx_delta(S_0, r, sigma, Del_TT, KK, delta_c_E, delta_p_E)
    call approx_gamma(S_0, r, sigma, Del_TT, KK, gamma_c_E, gamma_p_E)
    call approx_vega(S_0, r, sigma, Del_TT, KK, vega_c_E, vega_p_E)
    call approx_theta(S_0, r, sigma, Del_TT, KK, theta_c_E, theta_p_E)
    call approx_rho(S_0, r, sigma, Del_TT, KK, rho_c_E, rho_p_E)

    ! write output
    write(*,'(a)')'SENSITIVITY ANALYSIS FOR EUROPEAN OPTIONS'
    write(*,'(a)')'              Price       Delta       Gamma        Vega       Theta         Rho'
    write(*,'(a, 6f12.5)')'Call = ', c_E, delta_c_E, gamma_c_E, vega_c_E, theta_c_E, rho_c_E
    write(*,'(a, 6f12.5)')'Put  = ', p_E, delta_p_E, gamma_p_E, vega_p_E, theta_p_E, rho_p_E

    ! set up plot
    do in = 0, nplot

        S_plot(in) = 20d0 + dble(in)/10d0

        call black_scholes(S_plot(in), r, sigma, Del_TT, KK, p_plot(1, in), p_plot(2, in))
        call approx_delta(S_plot(in), r, sigma, Del_TT, KK, delta_plot(1, in), delta_plot(2, in))
        call approx_gamma(S_plot(in), r, sigma, Del_TT, KK, gamma_plot(1, in), gamma_plot(2, in))
        call approx_vega(S_plot(in), r, sigma, Del_TT, KK, vega_plot(1, in), vega_plot(2, in))
        call approx_theta(S_plot(in), r, sigma, Del_TT, KK, theta_plot(1, in), theta_plot(2, in))
        call approx_rho(S_plot(in), r, sigma, Del_TT, KK, rho_plot(1, in), rho_plot(2, in))

    enddo

    ! execute plot
    call plot(S_plot, p_plot(1, :), legend='European Call')
    call plot(S_plot, p_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Price')

    call plot(S_plot, delta_plot(1, :), legend='European Call')
    call plot(S_plot, delta_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Delta')

    call plot(S_plot, gamma_plot(1, :), legend='European Call')
    call plot(S_plot, gamma_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Gamma')

    call plot(S_plot, Vega_plot(1, :), legend='European Call')
    call plot(S_plot, Vega_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Vega (per %)')

    call plot(S_plot, theta_plot(1, :), legend='European Call')
    call plot(S_plot, theta_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Theta (per day)')

    call plot(S_plot, rho_plot(1, :), legend='European Call')
    call plot(S_plot, rho_plot(2, :), legend='European Put')
    call execplot(xlabel='Stock Price', ylabel='Rho (per %)')

end program
