!##############################################################################
! PROGRAM Bermuda
!
! ## Option pricing using binomial tree for Bermuda options
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
program Bermuda

    use toolbox
    implicit none
    real*8, parameter :: Del_TT_3 = 62d0/250d0   ! exercise date (in annualized values)
    real*8, parameter :: S_0 = 24d0              ! initial stock price
    real*8, parameter :: r = 0.04d0              ! annual interest rate
    real*8, parameter :: sigma = 0.10d0          ! standard deviation of stock returns
    real*8, parameter :: KK = 25d0               ! strike price
    integer, parameter :: TT_3 = 62              ! number of sub-periods (total)
    integer, parameter :: TT_2 = 40              ! number of sub-periods (second date)
    integer, parameter :: TT_1 = 20              ! number of sub-periods (first date)
    real*8 :: S_tn(0:TT_3, 0:TT_3), c_B(0:TT_3, 0:TT_3), p_B(0:TT_3, 0:TT_3)
    real*8 :: b, u, d, p, del_t
    integer :: it, n

    ! get length of sub-periods
    del_t = Del_TT_3/dble(TT_3)

    ! initialize parameters for binominal tree
    b = 0.5d0*(exp(-r*del_t) + exp((r+sigma**2)*del_t))
    u = b + sqrt(b**2-1d0)
    d = b - sqrt(b**2-1d0)
    p = (exp(r*del_t)-d)/(u-d)

    ! simulate stock price
    S_tn(0, 0) = S_0
    do it = 1, TT_3
        do n = 0, it
            S_tn(it, n) = S_0*u**n*d**(it-n)
        enddo
    enddo

    ! compute call and put price at maturity
    do n = 0, TT_3
        c_B(TT_3, n) = max(S_tn(TT_3, n) - KK, 0d0)
        p_B(TT_3, n) = max(KK - S_tn(TT_3, n), 0d0)
    enddo

    ! work recursively through the binomial tree
    do it = TT_3 - 1, 1, -1
        do n = 0, it

            ! check whether to exercise the option prior to maturity
            if(it == TT_1 .or. it == TT_2)then
                c_B(it, n) = max(S_tn(it, n) - KK, exp(-r*del_t) &
                                 *(p*c_B(it+1, n+1) + (1d0-p)*c_B(it+1, n)))
                p_B(it, n) = max(KK - S_tn(it, n), exp(-r*del_t) &
                                 *(p*p_B(it+1, n+1) + (1d0-p)*p_B(it+1, n)))
            else
                c_B(it, n) = exp(-r*del_t)*(p*c_B(it+1, n+1) &
                                            + (1d0-p)*c_B(it+1, n))
                p_B(it, n) = exp(-r*del_t)*(p*p_B(it+1, n+1) &
                                            + (1d0-p)*p_B(it+1, n))
            endif
        enddo
    enddo

    ! fair value of a bermuda option at the day of issue
    c_B(0, 0) = exp(-r*del_t)*(p*c_B(1, 1) + (1d0-p)*c_B(1, 0))
    p_B(0, 0) = exp(-r*del_t)*(p*p_B(1, 1) + (1d0 - p)*p_B(1, 0))

    write(*,'(a)')'BERMUDA OPTION PRICING USING BINOMIAL TREE'
    write(*,'(a, f12.5)')'Call Price = ', c_B(0, 0)
    write(*,'(a, f12.5)')'Put  Price = ', p_B(0, 0)

end program
