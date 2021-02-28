!##############################################################################
! PROGRAM StationaryDistribution
!
! ## Computes Stationary distribution of stochastic process.
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
include "sol_prog09_06m.f90"

program StationaryDistribution

    use globals
    use toolbox

    implicit none

    ! discretize stochastic process
    call discretize_AR(rho, 0d0, sigma, eta, pi)

    ! initialize distribution
    phi = 0d0
    phi(floor(dble(NS)/2d0)+1) = 1d0

    ! iterate until the distribution function converges
    do iter = 1, itermax

        phi_new = 0d0
        do is = 1, NS
            do is_p = 1, NS
                phi_new(is_p) = phi_new(is_p) + pi(is, is_p)*phi(is)
            enddo
        enddo

        ! get convergence level
        con_lev = maxval(abs(phi_new(:) - phi(:))/max(abs(phi(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            exit
        endif

        ! update distribution
        phi = phi_new
    enddo

    call output()

contains


    ! for creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: mu_u, sigma_u, rho_u, eta_l, eta_u
        real*8 :: eta_plot(0:nplot), cdf_plot(0:nplot)

        ! calculate unconditional distribution
        mu_u = sum(phi*eta)
        sigma_u = sum(phi*(eta-mu_u)**2)

        ! calculate autocorrelation
        rho_u = 0d0
        do is = 1, NS
            do is_p = 1, NS
                rho_u = rho_u + phi(is)*(eta(is)-mu_u)*(eta(is_p)-mu_u)*pi(is, is_p)
            enddo
        enddo
        rho_u = rho_u/sigma_u

        write(*,'(/a)')'        mu_u   sigma^2_u       rho_u       convergence'
        write(*,'(3f12.5, f18.8)')mu_u, sigma_u, rho_u, con_lev

        ! calculate cumulative distribution
        phi_cdf(1) = phi(1)/2d0
        do is = 2, NS
            phi_cdf(is) = phi_cdf(is-1) + (phi(is)+phi(is-1))/2d0
        enddo

        ! initialize lower and upper bound for plots
        eta_l = -3d0*sigma/(1d0-rho**2d0)
        eta_u = 3d0*sigma/(1d0-rho**2d0)

        ! calculate plot data for analytical cdf
        do is = 0, nplot
            eta_plot(is) = eta_l + (eta_u-eta_l)*dble(is)/dble(nplot)
            cdf_plot(is) = normalcdf(eta_plot(is), 0d0, sigma/(1d0-rho**2d0))
        enddo

        ! plot numerical and analytical cdf
        call plot(eta, phi_cdf, legend='numerical CDF (NS=5)')
        call plot(eta, phi_cdf, legend='numerical CDF (NS=15)')
        call plot(eta, phi_cdf, legend='numerical CDF (NS=25)')
        call plot(eta_plot, cdf_plot, legend='analytical CDF')
        call execplot(xlabel='eta', ylabel='Cumulative Density Function')

    end subroutine

end program
