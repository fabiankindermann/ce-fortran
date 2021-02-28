!##############################################################################
! PROGRAM household3
!
! ## The life cycle model with wage and interest rate risk
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
include "prog05_03m.f90"

program household3

    use globals
    use toolbox

    implicit none
    real*8 :: x(1+n_w*n_R)
    real*8 :: fret, low(1+n_w*n_R), up(1+n_w*n_R)
    integer :: j, iw, ir2, ir3, ic

    ! discretize w
    call log_normal_discrete(w, weight_w, mu_w, sig_w)

    ! discretize R
    call log_normal_discrete(R, weight_R, mu_R, sig_R)

    ! lower and upper border and initial guess
    low = 0d0
    up(1) = mu_w
    ic = 2
    do iw = 1, n_w
        do ir2 = 1, n_R
            up(ic) = R(ir2)*mu_w+w(iw)
            ic = ic+1
        enddo
    enddo
    x = up/2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! set up data for output
    do iw = 1, n_w
        do ir2 = 1, n_R
            do ir3 = 1, n_R
                wag(1,iw,ir2,ir3) = mu_w
                inc(1,iw,ir2,ir3) = mu_w
                sav(1,iw,ir2,ir3) = a(2,1,1)

                wag(2,iw,ir2,ir3) = w(iw)
                inc(2,iw,ir2,ir3) = w(iw)+R(ir2)*a(2,1,1)
                sav(2,iw,ir2,ir3) = a(3,iw,ir2)

                wag(3,iw,ir2,ir3) = 0d0
                inc(3,iw,ir2,ir3) = R(ir3)*a(3,iw,ir2)
                sav(3,iw,ir2,ir3) = 0d0
            enddo
        enddo
    enddo

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC    SAV'
    do j = 1, 3
        write(*,'(i4,4f7.2,a)')j,E(c(j,:,:,:)),E(wag(j,:,:,:)),&
            E(inc(j,:,:,:)),E(sav(j,:,:,:)),' (MEAN)'
        write(*,'(4x,4f7.2,a/)')Std(c(j,:,:,:)),Std(wag(j,:,:,:)),&
            Std(inc(j,:,:,:)),Std(sav(j,:,:,:)),' (STD)'
    enddo

    write(*,'(/2(a,f6.2))')' E(w) = ',sum(weight_w*w),'   Var(w) = ', &
        sum(weight_w*w**2)-sum(weight_w*w)**2
    write(*,'(2(a,f6.2))')' E(R) = ',sum(weight_R*R),'   Var(R) = ', &
        sum(weight_R*R**2)-sum(weight_R*R)**2

end program
