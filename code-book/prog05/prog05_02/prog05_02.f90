!##############################################################################
! PROGRAM household2
!
! ## The life cycle model with wage risk
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
include "prog05_02m.f90"

program household2

    use globals
    use toolbox

    implicit none
    real*8 :: x(1+n_w)
    real*8 :: fret, low(1+n_w), up(1+n_w)
    integer :: j, iw

    ! discretize w
    call log_normal_discrete(w, weight_w, mu_w, sig_w)

    ! lower and upper border and initial guess
    low = 0d0
    up(1) = mu_w
    up(2:1+n_w) = R*mu_w+w
    x = up/2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! set up data for output
    do iw = 1, n_w
        wag(1,iw) = mu_w
        inc(1,iw) = mu_w
        sav(1,iw) = a(2,1)

        wag(2,iw) = w(iw)
        inc(2,iw) = w(iw)+R*a(2,1)
        sav(2,iw) = a(3,iw)

        wag(3,iw) = 0d0
        inc(3,iw) = R*a(3,iw)
        sav(3,iw) = 0d0
    enddo

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC    SAV'
    do j = 1, 3
        write(*,'(i4,4f7.2,a)')j,E(c(j,:)),E(wag(j,:)),&
            E(inc(j,:)),E(sav(j,:)),' (MEAN)'
        write(*,'(4x,4f7.2,a/)')Std(c(j,:)),Std(wag(j,:)),&
            Std(inc(j,:)),Std(sav(j,:)),' (STD)'
    enddo

    write(*,'(/2(a,f6.2))')' E(w) = ',sum(weight_w*w),'   Var(w) = ', &
        sum(weight_w*w**2)-sum(weight_w*w)**2

end program
