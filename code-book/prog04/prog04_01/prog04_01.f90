!##############################################################################
! PROGRAM portfoliochoice
!
! ## Optimal portfolio choice, tangent portfolio, optimal risk exposure
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
program portfoliochoice

    use toolbox

    implicit none
    integer, parameter :: TT = 5        ! number of time periods
    integer, parameter :: NN = 3        ! number of stocks
    real*8, parameter :: gamma = 10d0   ! risk aversion
    real*8, parameter :: r_f = 0.05     ! the risk free rate

    real*8 :: omega_mv(NN), omega_tan(NN), omega(NN), omega_f
    real*8 :: mu(NN), sig(NN,NN), siginv(NN,NN), ID(NN)
    real*8 :: mu_mv, sig_mv, mu_tan, sig_tan, mu_p, sig_p

    integer :: ii
    real*8 :: a, b, c, x1(0:100), y1(0:100)

    ! get returns and their statistical properties
    call returns

    ! invert variance-covariance matrix
    siginv = lu_invert(sig)

    ! derive parameters
    a = dot_product(matmul(ID, siginv), ID)
    b = dot_product(matmul(ID, siginv), mu)
    c = dot_product(matmul(mu, siginv), mu)

    ! calculate minimum variance portfolio
    omega_mv = matmul(siginv,ID)/dot_product(matmul(ID,siginv),ID)

    ! check riskless interest rate
    if(dot_product(mu, omega_mv) < r_f)then
        write(*,'(a)') 'Riskless interest rate is too high!'
        stop
    endif

    ! compute tangent portfolio
    omega_tan = matmul(siginv, mu - r_f)

    ! get investors choice
    omega = omega_tan/gamma
    omega_f = 1d0 - sum(omega)

    ! normalize tangent portfolio shares to 1
    omega_tan = omega_tan/sum(omega_tan)

    ! compute returns and variances
    mu_mv   = dot_product(mu, omega_mv)
    sig_mv  = dot_product(omega_mv, matmul(sig, omega_mv))
    mu_tan  = dot_product(mu, omega_tan)
    sig_tan = dot_product(omega_tan, matmul(sig, omega_tan))
    mu_p    = R_f*omega_f + dot_product(mu, omega)
    sig_p   = dot_product(omega, matmul(sig, omega))

    ! print results
    write(*,'(a,3f8.3)')'Expected return risky asset',(mu(ii), ii = 1, 3)
    write(*,*)
    write(*,'(a,3f8.3)')'Variance-covariance matrix ',(sig(1,ii), ii = 1, 3)
    write(*,'(a,3f8.3)')'                           ',(sig(2,ii), ii = 1, 3)
    write(*,'(a,3f8.3)')'                           ',(sig(3,ii), ii = 1, 3)
    write(*,*)
    write(*,'(a)') 'Portfolio shares:   Riskless   Stock 1   Stock 2   Stock 3      mu_p     sig_p'
    write(*,'(a,10x,5f10.4)')'MV Portfolio      ', omega_mv(1), omega_mv(2), omega_mv(3), mu_mv, sig_mv
    write(*,'(a,10x,5f10.4)')'Tangent Portfolio ', omega_tan(1), omega_tan(2), omega_tan(3), mu_tan, sig_tan
    write(*,'(a,6f10.4)')'Optimal Portfolio ', omega_f, omega(1), omega(2), omega(3), mu_p, sig_p

    ! create plot data for efficient portfolios
    do ii = 0, 100
        x1(ii) = 1d0/1000d0*dble(ii+20)
    enddo
    y1 = (a*x1**2 - 2*b*x1+c)/(a*c - b**2)

    call plot(y1, x1)
    call plot((/1d0/a/), (/b/a/), marker=7, markersize=1.5d0, noline=.true., legend='MV Portfolio')
    call plot((/sig(1, 1)/), (/mu(1)/), marker=7, markersize=1.5d0, noline=.true., legend='Stock 1')
    call execplot(xlabel='Portfolio Variance', ylabel='Portfolio Mean Return')

contains


    ! get returns and statictical properties
    subroutine returns

        implicit none
        integer :: ii, ij
        real*8 :: odat(0:TT,NN), r(TT,NN)

        ! set vector of ones
        ID(:) = 1.0d0

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


