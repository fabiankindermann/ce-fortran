!##############################################################################
! PROGRAM household4
!
! ## The life cycle model with portfolio choice
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
include "prog05_04m.f90"

program household4

    use globals
    use toolbox

    implicit none
    real*8 :: x(2*(1+n_w*n_R))
    real*8 :: fret, low(2*(1+n_w*n_R)), up(2*(1+n_w*n_R))
    real*8 :: wag(3,n_w*n_R,n_R), inc(3,n_w*n_R,n_R)
    real*8 :: sav(3,n_w*n_R,n_R), alp(3,n_w*n_R,n_R)
    real*8 :: E_st(2), Var_st(2), rho_st
    integer :: j, iw, ir2, ir3, ic, iwR

    ! set tolerance level of minimization routine
    call settol_min(1d-14)

    ! discretize w and R in period 2
    call log_normal_discrete((/n_w, n_R/), wR, weight_wR, &
        (/mu_w, mu_R/), (/sig_w, sig_R/), rho_wR)

    ! discretize log(R) in period 3
    call log_normal_discrete(R, weight_R, mu_R, sig_R)

    ! lower and upper border and initial guess
    low = 0d0
    up(1) = mu_w
    up(2) = 1d0
    ic = 3
    iwR = 1
    do iw = 1, n_w
        do ir2 = 1, n_R
            up(ic) = wR(iwR,2)*mu_w+wR(iwR,1)
            up(ic+1) = 1d0
            ic = ic+2
            iwR = iwR+1
        enddo
    enddo
    x = up/2d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! set up data for output
    iwR = 1
    do iw = 1, n_w
        do ir2 = 1, n_R
            do ir3 = 1, n_R
                wag(1,iwR,ir3) = mu_w
                inc(1,iwR,ir3) = mu_w
                sav(1,iwR,ir3) = a(2,1)
                alp(1,iwR,ir3) = omega(1,1)*100d0

                wag(2,iwR,ir3) = wR(iwR,1)
                inc(2,iwR,ir3) = wR(iwR,1)+(Rf+omega(1,1)*(wR(iwR,2)-Rf))*a(2,1)
                sav(2,iwR,ir3) = a(3,iwR)
                alp(2,iwR,ir3) = omega(2,iwR)*100d0

                wag(3,iwR,ir3) = 0d0
                inc(3,iwR,ir3) = (Rf+omega(2,iwR)*(R(ir3)-Rf))*a(3,iwR)
                sav(3,iwR,ir3) = 0d0
                alp(3,iwR,ir3) = 0d0
            enddo
            iwR = iwR+1
        enddo
    enddo

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC    SAV  SHARE'
    do j = 1, 3
        write(*,'(i4,5f7.2,a)')j,E(c(j,:,:)),E(wag(j,:,:)),&
            E(inc(j,:,:)),E(sav(j,:,:)),E(alp(j,:,:)),' (MEAN)'
        write(*,'(4x,5f7.2,a/)')Std(c(j,:,:)),Std(wag(j,:,:)),&
            Std(inc(j,:,:))/E(inc(j,:,:)),Std(sav(j,:,:)),Std(alp(j,:,:)),' (STD)'
    enddo

    ! calculate w and R statistics
    E_st(1) = sum(wR(:,1)*weight_wR)
    E_st(2) = sum(wR(:,2)*weight_wR)
    Var_st(1) = sum(wR(:,1)**2*weight_wR)-E_st(1)**2
    Var_st(2) = sum(wR(:,2)**2*weight_wR)-E_st(2)**2
    rho_st = sum((wR(:,1)-E_st(1))*(wR(:,2)-E_st(2))*weight_wR)/&
        sqrt(max(Var_st(1)*Var_st(2), 1d-10))

    write(*,'(/2(a,f6.2))')' E(w)  = ',E_st(1),'   Var(w)  = ', &
        Var_st(1)
    write(*,'(3(a,f6.2))')' E(R2) = ',E_st(2),'   Var(R2) = ', &
        Var_st(2),'   rho(w, R2) = ',rho_st
    write(*,'(2(a,f6.2))')' E(R3) = ',sum(weight_R*R),'   Var(R3) = ', &
        sum(weight_R*R**2)-sum(weight_R*R)**2

end program


