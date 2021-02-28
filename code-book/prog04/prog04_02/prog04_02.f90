!##############################################################################
! PROGRAM portfoliochoice_short
!
! ## Optimal portfolio choice with short-selling constraints
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
include "prog04_02m.f90"

program portfoliochoice_short

    use globals
    use toolbox

    implicit none
    real*8 :: omega(NN), omega_f, mu_p, sig_p
    real*8 :: x_in(2*NN+1)
    logical :: check

    ! get returns and their statistical properties
    call returns

    ! initial guess for portfolio shares
    x_in(1:NN) = 0.2d0

    ! initial guess for multipliers
    x_in(NN+1:2*NN+1) = 0d0

    ! solve first oder conditions
    call fzero(x_in, focs, check)

    ! check for convergence
    if(check) stop 'Error: fzero did not converge'

    ! copy decision
    omega   = x_in(1:NN)
    omega_f = 1d0 - sum(omega)

    ! compute return and variance of portfolio
    mu_p   = R_f*omega_f + dot_product(mu, omega)
    sig_p  = dot_product(omega, matmul(sig, omega))

    ! print results
    write(*,*)
    write(*,'(a)') 'Portfolio shares:   Riskless   Stock 1   Stock 2   Stock 3      mu_p     sig_p'
    write(*,'(a,6f10.4)')'Optimal Portfolio ', omega_f, omega(1), omega(2), omega(3), mu_p, sig_p

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

    end subroutine

end program


