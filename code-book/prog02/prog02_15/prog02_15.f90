!##############################################################################
! PROGRAM distributions
!
! ## Probability density functions and cumulative distributions
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
program distributions

    use toolbox

    implicit none
    integer, parameter :: NN = 1000
    real*8 :: z(0:NN), dens(0:NN), dist(0:NN)
    real*8 :: a, b, mu, sigma, alpha, beta, p, q
    integer :: ii, n

    ! uniform distribution
    a = -1d0
    b = 1d0
    call grid_Cons_Equi(z, a-0.5d0, b+0.5d0)
    do ii = 0, NN
        dens(ii) = uniformPDF(z(ii), a, b)
        dist(ii) = uniformCDF(z(ii), a, b)
    enddo

    call plot(z, dens, legend="Density")
    call plot(z, dist, legend="Distribution")
    call execplot(title="Uniform Distribution")

    ! normal distribution
    mu = 1d0
    sigma = 0.25d0

    call grid_Cons_Equi(z, mu-4d0*sqrt(sigma), mu+4d0*sqrt(sigma))
    do ii = 0, NN
        dens(ii) = normalPDF(z(ii), mu, sigma)
        dist(ii) = normalCDF(z(ii), mu, sigma)
    enddo

    call plot(z, dens, legend="Density")
    call plot(z, dist, legend="Distribution")
    call execplot(title="Normal Distribution")

    ! log-normal distribution
    mu = 1d0
    sigma = 0.25d0

    call grid_Cons_Equi(z, 0d0, mu+5d0*sqrt(sigma))
    do ii = 0, NN
        dens(ii) = log_normalPDF(z(ii), mu, sigma)
        dist(ii) = log_normalCDF(z(ii), mu, sigma)
    enddo

    call plot(z, dens, legend="Density")
    call plot(z, dist, legend="Distribution")
    call execplot(title="Log-Normal Distribution")

    ! Gamma distribution
    alpha = 2d0
    beta = 1.25d0

    call grid_Cons_Equi(z, 0d0, 5d0)
    do ii = 0, NN
        dens(ii) = GammaPDF(z(ii), alpha, beta)
        dist(ii) = GammaCDF(z(ii), alpha, beta)
    enddo

    call plot(z, dens, legend="Density")
    call plot(z, dist, legend="Distribution")
    call execplot(title="Gamma Distribution")

    ! beta distribution
    p = 2d0
    q = 5d0

    call grid_Cons_Equi(z, -0.2d0, 1.2d0)
    do ii = 0, NN
        dens(ii) = betaPDF(z(ii), p, q)
        dist(ii) = betaCDF(z(ii), p, q)
    enddo

    call plot(z, dens, legend="Density")
    call plot(z, dist, legend="Distribution")
    call execplot(title="beta Distribution")

    ! bernoulli distribution
    p = 0.25d0
    do ii = 0, 1
        z(ii) = dble(ii)
        dens(ii) = bernoulliPDF(ii, p)
        dist(ii) = bernoulliCDF(ii, p)
    enddo

    call plot(z(0:1), dens(0:1), legend="Probabilities", marker=3, markersize=1.5d0, noline=.true.)
    call plot(z(0:1), dist(0:1), legend="Distribution", marker=3, markersize=1.5d0, noline=.true.)
    call execplot(title="Bernoulli Distribution")

    ! binomial distribution
    p = 0.25d0
    n = 12
    do ii = 0, n
        z(ii) = dble(ii)
        dens(ii) = binomialPDF(ii, n, p)
        dist(ii) = binomialCDF(ii, n, p)
    enddo

    call plot(z(0:n), dens(0:n), legend="Probabilities", marker=3, markersize=1.5d0, noline=.true.)
    call plot(z(0:n), dist(0:n), legend="Distribution", marker=3, markersize=1.5d0, noline=.true.)
    call execplot(title="Binomial Distribution")

end program
