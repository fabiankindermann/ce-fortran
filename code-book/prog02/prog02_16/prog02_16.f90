!##############################################################################
! PROGRAM simulation
!
! ## Simulation of realizations of random variables
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
program simulation

    use toolbox

    implicit none
    integer, parameter :: NN = 20000
    real*8 :: x(0:NN)
    real*8 :: a, b, mu, sigma, alpha, beta, p, q
    integer :: n

    ! uniform distribution
    a = -1d0
    b = 1d0
    call simulate_uniform(x, a, b)
    call plot_hist(x, 20)
    call execplot(title="Uniform Distribution")

    ! normal distribution
    mu = 1d0
    sigma = 0.25d0
    call simulate_normal(x, mu, sigma)
    call plot_hist(x, 20, left=mu-3d0*sqrt(sigma), right=mu+3d0*sqrt(sigma))
    call execplot(title="Normal Distribution")

    ! log-normal distribution
    mu = 1d0
    sigma = 0.25d0
    call simulate_log_normal(x, mu, sigma)
    call plot_hist(x, 20, left=0d0, right=mu+5d0*sqrt(sigma))
    call execplot(title="Log-Normal Distribution")

    ! Gamma distribution
    alpha = 2d0
    beta = 1.25d0
    call simulate_Gamma(x, alpha, beta)
    call plot_hist(x, 20)
    call execplot(title="Gamma Distribution")

    ! beta distribution
    p = 2d0
    q = 5d0
    call simulate_beta(x, p, q)
    call plot_hist(x, 20)
    call execplot(title="beta Distribution")

    ! bernoulli distribution
    p = 0.25d0
    call simulate_bernoulli(x, p)
    call plot_hist(x, 2)
    call execplot(title="Bernoulli Distribution")

    ! binomial distribution
    p = 0.25d0
    n = 12
    call simulate_binomial(x, n, p)
    call plot_hist(x, n+1, left=-0.5d0, right=dble(n)+0.5d0)
    call execplot(title="Binomial Distribution")

end program
