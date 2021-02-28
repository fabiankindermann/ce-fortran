!##############################################################################
! PROGRAM portfoliochoice_MCM
!
! ## Portfolio choice with Monte Carlo minimization (MCM)
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
program portfoliochoice_MCM

    use toolbox

    implicit none
    integer, parameter :: TT = 5        ! number of time periods
    integer, parameter :: NN = 3        ! number of stocks
    integer, parameter :: NP = 100      ! number of mc portfolios
    real*8, parameter :: gamma = 10d0   ! risk aversion
    real*8, parameter :: r_f = 0.05     ! the risk free rate

    real*8 :: omega(0:NN)
    real*8 :: mu(0:NN), sig(0:NN,0:NN)
    real*8 :: mu_p, sig_p

    real*8 :: omega_mc(0:NN,NP), mu_mc(NP), sig_mc(NP)
    integer :: ii, ik, maxi

    ! get returns and their statistical properties
    call returns

    ! simulate initial portfolios
    do ii = 1, NP
        call simulate_uniform(omega_mc(:, ii), 0d0, 1d0)
        omega_mc(:, ii) = omega_mc(:, ii)/sum(omega_mc(:, ii))
    enddo

    ! perform updating steps
    do ik = 1, 1000

        ! calculate means and returns of the portfolios
        do ii = 1, NP
            mu_mc(ii)  = dot_product(mu, omega_mc(:, ii))
            sig_mc(ii) = dot_product(omega_mc(:, ii), matmul(sig, omega_mc(:, ii)))
        enddo

        ! get portfolio that maximizes investors utility
        maxi = maxloc(mu_mc - gamma/2d0*sig_mc, 1)

        ! store maximum
        omega = omega_mc(:, maxi)

        ! simulate an updated set of portfolios around maximum
        do ii = 1, NP-1

            ! get new normalized portfolio
            call simulate_uniform(omega_mc(:, ii), 0d0, 1d0)
            omega_mc(:, ii) = omega_mc(:, ii)/sum(omega_mc(:, ii))

            ! create linear combination between old and new
            omega_mc(:, ii) = 0.9d0*omega + 0.1d0*omega_mc(:, ii)
        enddo

        ! store maximum in last entry
        omega_mc(:, NP) = omega
    enddo

    ! calculate return and variance
    mu_p   = dot_product(mu, omega)
    sig_p  = dot_product(omega, matmul(sig, omega))

    ! print results
    write(*,'(a)') 'Portfolio shares:   Riskless   Stock 1   Stock 2   Stock 3      mu_p     sig_p'
    write(*,'(a,6f10.4)')'Optimal Portfolio ', omega(0), omega(1), omega(2), omega(3), mu_p, sig_p

contains


    ! get returns and statictical properties
    subroutine returns

        implicit none
        integer :: ii, ij
        real*8 :: odat(0:TT,NN), r(TT,NN)

        ! initialize price data
        odat(0,:) = (/ 1.00d0, 2.00d0, 3.00d0 /)
        odat(1,:) = (/ 1.02d0, 2.65d0, 2.80d0 /)
        odat(2,:) = (/ 1.17d0, 2.40d0, 4.50d0 /)
        odat(3,:) = (/ 1.08d0, 2.70d0, 4.20d0 /)
        odat(4,:) = (/ 1.16d0, 2.85d0, 3.20d0 /)
        odat(5,:) = (/ 1.26d0, 2.75d0, 4.20d0 /)

        ! compute return per period
        do ii = 1, NN
            do ij = 1, TT
                r(ij, ii) = odat(ij, ii)/odat(ij-1, ii)-1d0
            enddo
        enddo

        ! compute expected value
        do ii = 1, NN
            mu(ii) = sum(r(:, ii))/dble(TT)
        enddo

        ! compute covariance matrix
        do ii = 1, NN
            do ij = 1, NN
                sig(ii, ij) = dot_product(r(:, ii)-mu(ii), r(:, ij)-mu(ij))
            enddo
        enddo
        sig = sig/dble(TT)

        ! store risk free rate and variance
        mu(0) = r_f
        sig(0, :) = 0d0
        sig(:, 0) = 0d0

    end subroutine

end program

