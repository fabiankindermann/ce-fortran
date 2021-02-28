!##############################################################################
! PROGRAM Option_Pricing
!
! ## Option pricing using binomial tree and the Black-Scholes formula
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
program Option_Pricing

    use toolbox

    implicit none
    real*8, parameter :: Del_TT = 62d0/250d0     ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 25d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    integer, parameter :: TT = 200               ! number of sub-periods

    real*8 :: S(0:TT), z(TT), t(0:TT)
    real*8 :: del_t, d_1, d_2, b, u, d, p
    real*8 :: c_E, p_E, c_A(0:TT,0:TT), p_A(0:TT,0:TT)
    real*8 :: S_tn(0:TT, 0:TT)
    integer :: it, n

    ! get length of sub-periods
    del_t = Del_TT/dble(TT)

    ! generate standard normal random numbers
    call simulate_normal(z, 0d0, 1d0)

    ! simulate stock price movement
    S(0) = S_0
    do it = 1, TT
        t(it) = dble(it)*del_t*250d0
        S(it) = S(it-1)*exp((r-sigma**2/2d0)*del_t + sigma*sqrt(del_t)*z(it))
    enddo

    ! plot simulated stock price
    call plot(t, S)
    call execplot()

    ! use Black-Scholes formula for call and put price
    d_1 = (log(S_0/KK) + (r+sigma**2/2)*Del_TT)/(sigma*sqrt(Del_TT))
    d_2 = d_1 - sigma*sqrt(Del_TT)
    c_E =  S_0*normalCDF( d_1, 0d0, 1d0) - KK*exp(-r*Del_TT)*normalCDF( d_2, 0d0, 1d0)
    p_E = -S_0*normalCDF(-d_1, 0d0, 1d0) + KK*exp(-r*Del_TT)*normalCDF(-d_2, 0d0, 1d0)

    write(*,'(a)')'OPTION PRICING USING BLACK-SCHOLES FORMULA'
    write(*,'(a, f12.5)')'Call Price = ', c_E
    write(*,'(a, f12.5)')'Put  Price = ', p_E
    write(*,*)

    ! initialize parameters for a binominal tree
    b = 0.5d0*(exp(-r*del_t) + exp((r + sigma**2)*del_t))
    u = b + sqrt(b**2 - 1d0)
    d = b - sqrt(b**2 - 1d0)
    p = (exp(r*del_t) - d)/(u - d)

    ! determine stock price realizations
    S_tn(0, 0) = S_0
    do it = 1, TT
        do n = 0, it
            S_tn(it, n) = S_0*u**n*d**(it-n)
        enddo
    enddo

    ! compute call and put price of a European option
    c_E = 0d0
    p_E = 0d0
    do n = 0, TT
        c_E = c_E + binomialPDF(n, TT, p)*max(S_tn(TT, n) - KK,0d0)
        p_E = p_E + binomialPDF(n, TT, p)*max(KK - S_tn(TT, n),0d0)
    enddo
    c_E = c_E*exp(-r*Del_TT)
    p_E = p_E*exp(-r*Del_TT)

    ! compute call and put price of American option in exercise period
    do n = 0, TT
        c_A(TT, n) = max(S_tn(TT, n) - KK, 0d0)
        p_A(TT, n) = max(KK - S_tn(TT, n), 0d0)
    enddo

    ! recursive computation of premiums
    do it = TT-1, 1, -1
        do n = 0, it
            c_A(it, n) = max(S_tn(it, n) - KK, exp(-r*del_t)*(p*c_A(it+1, n+1) + (1d0-p)*c_A(it+1, n)))
            p_A(it, n) = max(KK - S_tn(it, n), exp(-r*del_t)*(p*p_A(it+1, n+1) + (1d0-p)*p_A(it+1, n)))
        enddo
    enddo
    c_A(0, 0) = exp(-r*del_t)*(p*c_A(1, 1) + (1d0-p)*c_A(1, 0))
    p_A(0, 0) = exp(-r*del_t)*(p*p_A(1, 1) + (1d0-p)*p_A(1, 0))

    write(*, '(a)')'OPTION PRICING USING BINOMIAL TREE'
    write(*, '(a)')'                 European    American'
    write(*, '(a, 2f12.5)')'Call Price = ', c_E, c_A(0, 0)
    write(*, '(a, 2f12.5)')'Put  Price = ', p_E, p_A(0, 0)

end program
