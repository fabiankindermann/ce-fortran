!##############################################################################
! PROGRAM bequest
!
! ## The life cycle model with bequests
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog05_06m.f90"

program bequest

    use globals
    use toolbox

    implicit none
    real*8 :: x(3*n_w+2)
    real*8 :: fret, low(3*n_w+2), up(3*n_w+2)
    integer :: j, iw, ic

    ! set tolerance level of minimization routine
    call settol_min(1d-14)

    ! discretize log(wage)
    call log_normal_discrete(w, weight_w, mu_w, sig_w)

    ! calculate annuity factors
    p_a(1) = psi(2)/R+psi(2)*psi(3)/R**2
    p_a(2) = psi(3)/R

    ! lower and upper border and initial guess
    low = 0d0
    up(1) = mu_w
    up(2) = 1d0
    ic = 3
    do iw = 1, n_w
        up(ic) = R*mu_w + w(iw)
        up(ic + 1) = 1d0
        up(ic + 2) = R*(R*mu_w + w(iw)) + pen
        ic = ic+3
    enddo
    x = 0.2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)
    call fminsearch(x, fret, low, up, utility)

    ! set up data for output
    do iw = 1, n_w
        wag(1, iw) = mu_w
        inc(1, iw) = mu_w
        sav(1, iw) = a(2, 1)*(1d0 - omega(1, 1))
        alp(1, iw) = a(2, 1)*omega(1, 1)

        wag(2, iw) = w(iw)
        inc(2, iw) = w(iw) + R*(1d0 - omega(1, 1))*a(2, 1)+ &
                    omega(1, 1)*a(2, 1)/p_a(1)
        sav(2, iw) = a(3, iw)*(1d0 - omega(2, iw))
        alp(2, iw) = a(3, iw)*omega(2, iw)

        wag(3, iw) = 0d0
        inc(3, iw) = R*(1d0 - omega(2, iw))*a(3, iw) + omega(2, iw)*a(3, iw)/p_a(2) + &
                    omega(1, 1)*a(2, 1)/p_a(1) + pen
        sav(3, iw) = a(4, iw)
        alp(3, iw) = 0d0
    enddo


    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC   SREG   SANN'
    do j = 1, 3
        write(*,'(i4,5f7.2,a)')j,E(c(j,:)),E(wag(j,:)),&
            E(inc(j,:)),E(sav(j,:)),E(alp(j,:)),' (MEAN)'
        write(*,'(4x,5f7.2,a/)')Std(c(j,:)),Std(wag(j,:)),&
            Std(inc(j,:)),Std(sav(j,:)),Std(alp(j,:)),' (STD)'
    enddo

    write(*,'(/2(a,f6.2))')' E(w) = ',sum(weight_w*w),'   Var(w) = ', &
        sum(weight_w*w**2)-sum(weight_w*w)**2

    write(*,'(a,f6.2)')' omega1 = ',omega(1,1)
    write(*,'(a,f6.2)')' omega2 = ',E(omega(2,:))

    write(*,'(/a,f7.3)')'util   = ',-fret

end program
