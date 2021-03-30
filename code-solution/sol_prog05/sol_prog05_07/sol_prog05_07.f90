!##############################################################################
! PROGRAM equityannuity
!
! ## The life cycle model with annuities and portfolio choice
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog05_07m.f90"

program equityannuity

    use globals
    use toolbox

    implicit none
    real*8 :: x(3*(1+n_w*n_R))
    real*8 :: fret, low(3*(1+n_w*n_R)), up(3*(1+n_w*n_R))
    integer :: j, iw, ir2, ir3, ic, iwR

    ! set tolerance level of minimization routine
    call settol_min(1d-10)

    ! discretize w and R in period 2
    call log_normal_discrete((/n_w, n_R/), wR, weight_wR, &
        (/mu_w, mu_R/), (/sig_w, sig_R/), rho_wR)

    ! discretize log(R) in period 3
    call log_normal_discrete(R, weight_R, mu_R, sig_R)

    ! calculate annuity factors
    pa(1) = (psi(2)/hatR + psi(2)*psi(3)/hatR**2) * (1d0 + xi)
    pa(2) = (psi(3)/hatR) * (1d0 + xi)

    ! lower and upper border and initial guess
    low = 0d0
    up(1) = mu_w
    up(2) = 1.0d0
    up(3) = 1.0d0
    ic = 4
    iwR = 1
    do iw = 1, n_w
        do ir2 = 1, n_R
            up(ic) = wR(iwR, 2)*mu_w+wR(iwR, 1)
            up(ic + 1) = 1.0d0
            up(ic + 2) = 1.0d0
            ic = ic + 3
            iwR = iwR + 1
        enddo
    enddo
    x = up/10d0

    ! minimization routine
    call fminsearch(x, fret, low, up, utility)

    ! set up data for output
    iwR = 1
    do iw = 1, n_w
        do ir2 = 1, n_R
            do ir3 = 1, n_R
                wag(1, iwR, ir3) = mu_w
                inc(1, iwR, ir3) = mu_w
                sh(1) = (1d0 - omega_e(1, 1))*(1d0 - omega_a(1, 1))
                sh(2) = omega_e(1, 1)
                sh(3) = (1d0 - omega_e(1, 1))*omega_a(1, 1)
                sav_b(1, iwR, ir3) = a(2, 1)*sh(1)
                sav_e(1, iwR, ir3) = a(2, 1)*sh(2)
                sav_a(1, iwR, ir3) = a(2, 1)*sh(3)
                wag(2, iwR, ir3) = wR(iwR, 1)
                inc(2, iwR, ir3) = wR(iwR, 1)+(sh(1)*Rf + sh(2)*wR(iwR,2) + sh(3)/pa(1))*a(2, 1)

                sh(1) = (1d0 - omega_e(2, iwR))*(1d0 - omega_a(2, iwR))
                sh(2) = omega_e(2,iwR)
                sh(3) = (1d0 - omega_e(2, iwR))*omega_a(2, iwR)

                sav_b(2, iwR, ir3) = a(3, iwR)*sh(1)
                sav_e(2, iwR, ir3) = a(3, iwR)*sh(2)
                sav_a(2, iwR, ir3) = a(3, iwR)*sh(3)
                wag(3, iwR, ir3) = 0d0
                inc(3, iwR, ir3) = (sh(1)*Rf + sh(2)*R(ir3)+sh(3)/pa(2))*a(3, iwR) + &
                             (1d0-omega_e(1, 1))*omega_a(1, 1)*a(2, 1)/pa(1) + pen
                sav_b(3, iwR, ir3) = 0d0
                sav_e(3, iwR, ir3) = 0d0
                sav_a(3, iwR, ir3) = 0d0
            enddo
            iwR = iwR + 1
        enddo
    enddo

    ! output
    write(*,'(/a/)')' AGE   CONS   WAGE    INC   SBON   SEQU   SANN'
    do j = 1, 3
        write(*,'(i4,6f7.2,a)')j,E(c(j,:,:)),E(wag(j,:,:)),&
            E(inc(j,:,:)),E(sav_b(j,:,:)),E(sav_e(j,:,:)),E(sav_a(j,:,:)),' (MEAN)'
        write(*,'(4x,6f7.2,a/)')Std(c(j,:,:)),Std(wag(j,:,:)),&
            Std(inc(j,:,:)),Std(sav_b(j,:,:)),Std(sav_e(j,:,:)),Std(sav_a(j,:,:)),' (STD)'
    enddo

    ! calculate statistics for w and R
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

    write(*,'(/a,f6.3)')' util = ', -fret

end program
