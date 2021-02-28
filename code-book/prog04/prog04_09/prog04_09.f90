!##############################################################################
! PROGRAM mortality_riskmanagement
!
! ## Mortality risk management with Monte Carlo simulations
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
include "prog04_09m.f90"

program mortality_riskmanagement

    use toolbox
    use globals

    implicit none

    ! initialize mortality data
    call init_data


    !########################################################
    ! SIMULATE DIFFERENT MORTAILITY TABLES
    !########################################################

    ! MC simulations of mortality tables
    do ik = 1, KK

        ! simulate innovations to the random walk
        call simulate_normal(eps, 0d0, sigma_eps)

        ! calculate common risk-factor
        kappa(0) = 0d0
        do it = 1, TT
            kappa(it) = kappa(it-1) + c + eps(it)
        enddo

        ! determine conditional mortality rates (Lee-Carter model)
        do it = 0, TT
            q(:, it, ik) = exp(alpha + beta*kappa(it))
        enddo
    enddo

    ! store current probabilities in entry zero
    do it = 1, TT
        q(:, it, 0) = q(:, 0, 0)
    enddo

    ! calculate unconditional survival probabilities and life expectancy
    do ik = 0, KK
        do ix = 1, XX
            LE(ix, ik) = 0d0
            PP(ix, 0, ik) = 1d0
            do it = 1, XX - ix
                PP(ix, it, ik) = PP(ix, it-1, ik)*(1d0 - q(ix+it-1, it-1, ik))
                LE(ix, ik) = LE(ix, ik) + PP(ix, it, ik)*q(ix+it, it, ik)*dble(it)
            enddo
        enddo
    enddo

    ! determine path with highest and lowest life expectancy (at age 20-24)
    ik_min = minloc(LE(1, 1:KK), 1)
    ik_max = maxloc(LE(1, 1:KK), 1)

    ! determine plotting data
    do it = 0, TT
        age(it) = 20d0 + 5d0*dble(it)
        PP_av(it)  = sum(PP(1, it, 1:KK))/dble(KK)
        PP_min(it) = PP(1, it, ik_min)
        PP_max(it) = PP(1, it, ik_max)
    enddo
    age(TT+1)    = 100d0
    PP_av(TT+1)  = 0d0
    PP_min(TT+1) = 0d0
    PP_max(TT+1) = 0d0

    ! create plot of unconditional survival probabilities
    call plot(age, PP_av, color='blue', legend='Average')
    call plot(age, PP_min, color='blue', legend='Minimum', marker=6, noline=.true.)
    call plot(age, PP_max, color='blue', legend='Maximum', marker=4, noline=.true.)
    call execplot(xlabel='Age x', ylabel='Uncond. Survival Probability')

    ! print results
    write(*,*)
    write(*,'(a/)')'ESTIMATES FOR LIFE EXPECTANCY AT AGE 20'
    write(*,'(a,f10.3)')'Average = ', 20d0 + 5d0*sum(LE(1, :))/dble(KK)
    write(*,'(a,f10.3)')'Minimum = ', 20d0 + 5d0*LE(1, ik_min)
    write(*,'(a,f10.3)')'Maximum = ', 20d0 + 5d0*LE(1, ik_max)
    write(*,*)


    !########################################################
    ! CALCULATE PREMIA OF INSURANCE POLICIES
    !########################################################

    ! premium calculation for different insurance polices
    do ii = 1, NN
        ix = ins(ii)%ix
        im = ins(ii)%im
        id = ins(ii)%id

        if(ins(ii)%itype == 'AI')then
            p_m(ii) = AI(ix, im, id, 0)*(1d0 + xi_1)
        else
            p_m(ii) = LI(ix, im, id, 0)*(1d0 - xi_2)
        endif
    enddo

    ! get the total sum of premia payed for insuance contracts
    sum_premia = sum(ins(:)%n_ins*p_m*ins(:)%BB)

    ! determine shares in the insurer's portfolio
    omega = ins(:)%n_ins*p_m*ins(:)%BB/sum_premia

    ! print results
    write(*,'(/a/)')'PORTFOLIO OF THE INSURER'
    write(*,'(a)')'  #       x   m   d  Benefits     Price    Weight'
    do ii = 1, NN
        write(*, '(i3,a4,3i4,f10.1,2f10.3)')ii, ins(ii)%itype, &
                ins(ii)%ix, ins(ii)%im, ins(ii)%id, ins(ii)%BB, p_m(ii), omega(ii)*100d0
    enddo
    write(*, '(/a,f12.1/)')'Sum of premia payed (in million Euro)', &
        sum_premia/1000000d0


    !########################################################
    ! SIMULATE BENEFIT-PREMIUM RATIOS
    !########################################################

    ! calculate benefit-premium ratios for all simulations
    do ii = 1, NN
        ix = ins(ii)%ix
        im = ins(ii)%im
        id = ins(ii)%id

        if(ins(ii)%itype == 'AI')then
            do ik = 1, KK
                l(ii, ik) = AI(ix, im, id, ik)/p_m(ii)
            enddo
        else
            do ik = 1, KK
                l(ii, ik) = LI(ix, im, id, ik)/p_m(ii)
            enddo
        endif
    enddo

    ! sum up portfolio weighted benefit-premium ratios
    do ik = 1, KK
        l_p(ik) = sum(omega*l(1:NN, ik))
    enddo

    ! determine mean of benefit-premium ratios
    do ii = 1, NN
        mu(ii) = sum(l(ii,:))/dble(KK)
    enddo
    mu_p = sum(omega*mu)

    ! compute variance-covariance matrix of benefit-premium ratios
    sig_p = 0d0
    do ii = 1, NN
        do ij = 1, NN
            sig(ii, ij) = sum((l(ii,:)-mu(ii))*(l(ij,:)-mu(ij)))/dble(KK-1)
            sig_p = sig_p + sig(ii,ij)*omega(ii)*omega(ij)
        enddo
    enddo

    ! compute VaR and CVaR for portfolio
    call sort(l_p)
    ii_VaR = ceiling(0.95d0*KK)
    VaR    = l_p(ii_VaR)
    CVaR   = sum(l_p(ii_VaR:KK))/dble(KK-ii_VaR+1)

    ! print results
    write(*,'(/a/)')'PORTFOLIO CHARACTERISTICS'
    write(*,'(a, 5f10.3/)')'mu    = ', mu
    write(*,'(a, 5f10.3)')'sig   = ', sig(ii, :)
    do ii = 2, NN
        write(*,'(8x,5f10.3)')sig(ii, :)
    enddo
    write(*,'(/2(a,f10.3,4x))')'mu_p  = ', mu_p, 'sig_p = ', sig_p
    write(*,'(2(a,f10.3,4x)/)')'VaR   = ',  VaR, 'CVaR  = ',  CVar


    !########################################################
    ! COMPUTE OPTIMAL PORTFOLIO
    !########################################################

    ! initial guess for portfolio shares
    x_in(1:NN) = 1d0/dble(NN)

    ! initial guess for multipliers
    x_in(NN+1:2*NN+1) = 0d0

    ! solve first oder conditions
    call fzero(x_in, focs, check)

    ! check for convergence
    if(check) stop 'Error: fzero did not converge'

    ! copy decision
    omega   = x_in(1:NN)


    !########################################################
    ! GET PROPERTIES OF OPTIMAL PORTFOLIO
    !########################################################

    ! calculate new portfolio mean, variance and VaR's
    mu_p  = sum(omega*mu)
    sig_p = dot_product(matmul(omega, sig), omega)
    n_new = omega*sum_premia/(p_m*ins(:)%BB)

    ! sum up portfolio weighted benefit-premium ratios
    do ik = 1, KK
        l_p(ik) = sum(omega*l(1:NN, ik))
    enddo

    ! compute VaR and CVaR for portfolio
    call sort(l_p)
    ii_VaR = ceiling(0.95d0*KK)
    VaR    = l_p(ii_VaR)
    CVaR   = sum(l_p(ii_VaR:KK))/dble(KK-ii_VaR+1)

    ! print results
    write(*,'(/a/)')'OPTIMAL PORTFOLIO'
    write(*,'(a, 5f10.3)')'omega = ', omega
    write(*,'(a, 5f10.3)')'n     = ', n_new
    write(*,'(/2(a,f10.3,4x))')'mu_p  = ', mu_p, 'sig_p = ', sig_p
    write(*,'(2(a,f10.3,4x)/)')'VaR   = ',  VaR, 'CVaR  = ',  CVar

contains


    ! expected present value of annuity insurance benefits
    function AI(ix, im, id, ik)

        implicit none
        integer, intent(in) :: ix, im, id, ik
        real*8 :: AI
        integer :: it

        AI = 0d0
        do it = im, im+id-1
            AI = AI + PP(ix, it, ik)/(1+r)**it
        enddo

    end function AI


    ! expected present value of life insurance benefits
    function LI(ix, im, id, ik)

        implicit none
        integer, intent(in) :: ix, im, id, ik
        real*8 :: LI
        integer :: it

        LI = 0d0
        do it = im, im+id-1
            LI = LI + PP(ix, it, ik)*q(ix+it, it, ik)/(1+r)**(it+1)
        enddo

    end function LI


    ! initialize mortality and portfolio data
    subroutine init_data()

        implicit none

        q = 0d0

        ! mortality probabilities in the base period (t = 0)
        q(1,  0, :) = 0.00186d0
        q(2,  0, :) = 0.00222d0
        q(3,  0, :) = 0.00268d0
        q(4,  0, :) = 0.00389d0
        q(5,  0, :) = 0.00632d0
        q(6,  0, :) = 0.01123d0
        q(7,  0, :) = 0.01957d0
        q(8,  0, :) = 0.02987d0
        q(9,  0, :) = 0.04483d0
        q(10, 0, :) = 0.06528d0
        q(11, 0, :) = 0.10007d0
        q(12, 0, :) = 0.17166d0
        q(13, 0, :) = 0.28709d0
        q(14, 0, :) = 0.46296d0
        q(15, 0, :) = 0.66523d0
        q(16, :, :) = 1.00000d0

        ! parameters of the Lee-Carter Model !!! Here 0.02 in first entry
        alpha = log(q(:, 0, 1))
        beta  = (/ 0.02d0, 0.12d0, 0.12d0, 0.11d0, 0.10d0, 0.08d0, 0.07d0, 0.06d0, &
                   0.05d0, 0.05d0, 0.05d0, 0.05d0, 0.05d0, 0.04d0, 0.03d0, 0.00d0 /)

        ! define different insurance types
        ins(1) = insurance_type('AI', 9, 1, 7,  50000d0, 1000d0)
        ins(2) = insurance_type('AI',10, 0, 7,  50000d0, 1000d0)
        ins(3) = insurance_type('LI', 6, 0, 7, 500000d0, 1000d0)
        ins(4) = insurance_type('LI', 8, 0, 7, 500000d0, 1000d0)
        ins(5) = insurance_type('LI', 8, 0, 9, 500000d0, 1000d0)

    end subroutine init_data

end program
