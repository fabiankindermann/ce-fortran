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
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    real*8, parameter :: delta = 1d-5            ! approximation level for difference quotient
    integer, parameter :: nplot = 100            ! number of gridpoints for plot

contains

    subroutine black_scholes(S_in, r_in, sigma_in, Del_TT_in, KK_in, c_out, p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: c_out, p_out
        real*8 :: d_1, d_2

        ! compute call and put price of European option
        d_1 = (log(S_in/KK_in) + (r_in+sigma_in**2/2d0)*Del_TT_in)/sigma_in/sqrt(Del_TT_in)
        d_2 = d_1 - sigma_in*sqrt(Del_TT_in)
        c_out = S_in*normalCDF(d_1, 0d0, 1d0) - KK_in*exp(-r_in*Del_TT_in)*normalCDF(d_2, 0d0, 1d0)
        p_out = -S_in*normalCDF(-d_1, 0d0, 1d0) + KK_in*exp(-r_in*Del_TT_in)*normalCDF(-d_2, 0d0, 1d0)

    end subroutine

    subroutine approx_delta(S_in, r_in, sigma_in, Del_TT_in, KK_in, &
                            delta_c_out, delta_p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: delta_c_out, delta_p_out
        real*8 :: c_temp(2), p_temp(2)

        ! compute delta of an option
        call black_scholes(S_in, r_in, sigma_in, Del_TT_in, &
                           KK_in, c_temp(1), p_temp(1))
        call black_scholes(S_in + delta, r_in, sigma_in, Del_TT_in, &
                           KK_in, c_temp(2), p_temp(2))

        delta_c_out = (c_temp(2) - c_temp(1))/delta
        delta_p_out = (p_temp(2) - p_temp(1))/delta

    end subroutine

    subroutine approx_gamma(S_in, r_in, sigma_in, Del_TT_in, KK_in, &
                            gamma_c_out, gamma_p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: gamma_c_out, gamma_p_out
        real*8 :: delta_c_temp(2), delta_p_temp(2)

        ! compute delta of an option
        call approx_delta(S_in, r_in, sigma_in, Del_TT_in, KK_in, delta_c_temp(1), delta_p_temp(1))
        call approx_delta(S_in + delta, r_in, sigma_in, Del_TT_in, KK_in, delta_c_temp(2), delta_p_temp(2))

        gamma_c_out = (delta_c_temp(2) - delta_c_temp(1))/delta
        gamma_p_out = (delta_p_temp(2) - delta_p_temp(1))/delta

    end subroutine

    subroutine approx_vega(S_in, r_in, sigma_in, Del_TT_in, KK_in, &
                           vega_c_out, vega_p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: vega_c_out, vega_p_out
        real*8 :: c_temp(2), p_temp(2)

        ! compute vega of an option
        call black_scholes(S_in, r_in, sigma_in, Del_TT_in, KK_in, c_temp(1), p_temp(1))
        call black_scholes(S_in, r_in, sigma_in + delta, Del_TT_in, KK_in, c_temp(2), p_temp(2))

        vega_c_out = (c_temp(2) - c_temp(1))/delta/100d0
        vega_p_out = (p_temp(2) - p_temp(1))/delta/100d0

    end subroutine

    subroutine approx_theta(S_in, r_in, sigma_in, Del_TT_in, KK_in, &
                            theta_c_out, theta_p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: theta_c_out, theta_p_out
        real*8 :: c_temp(2), p_temp(2)

        ! compute theta of an option
        call black_scholes(S_in, r_in, sigma_in, Del_TT_in, KK_in, c_temp(1), p_temp(1))
        call black_scholes(S_in, r_in, sigma_in, Del_TT_in - delta, KK_in, c_temp(2), p_temp(2))

        theta_c_out = (c_temp(2) - c_temp(1))/delta/250d0
        theta_p_out = (p_temp(2) - p_temp(1))/delta/250d0

    end subroutine

    subroutine approx_rho(S_in, r_in, sigma_in, Del_TT_in, KK_in, &
                          rho_c_out, rho_p_out)

        implicit none

        real*8, intent(in) :: S_in, r_in, sigma_in, Del_TT_in, KK_in
        real*8, intent(out) :: rho_c_out, rho_p_out
        real*8 :: c_temp(2), p_temp(2)

        ! compute rho of an option
        call black_scholes(S_in, r_in, sigma_in, Del_TT_in, KK_in, c_temp(1), p_temp(1))
        call black_scholes(S_in, r_in + delta, sigma_in, Del_TT_in, KK_in, c_temp(2), p_temp(2))

        rho_c_out = (c_temp(2) - c_temp(1))/delta/100d0
        rho_p_out = (p_temp(2) - p_temp(1))/delta/100d0

    end subroutine

end module
