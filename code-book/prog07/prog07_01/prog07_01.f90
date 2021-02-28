!##############################################################################
! PROGRAM TRVL_OLG
!
! ## The OLG model with variable labor supply
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
program TRVL_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 24
    integer, parameter :: JJ      = 3
    integer, parameter :: JR      = 3
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  nu      = 1.5d0
    real*8, parameter ::  rho     = 0.6d0
    real*8, parameter ::  erho    = 1d0 - 1d0/rho
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.3d0
    integer, parameter :: itermax = 200

    ! model variables
    real*8 :: w(0:TT), r(0:TT), mu(JJ, 0:TT), wn(JJ, 0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT), pen(JJ,0:TT)
    real*8 :: by(0:TT), kappa(0:TT), n_p(0:TT), gy, tauk(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BA(0:TT), BF(0:TT), TB(0:TT)
    real*8 :: Tpen(0:TT), TXR(0:TT), h(JJ), m(JJ,0:TT)
    real*8 :: a(JJ,0:TT), c(JJ,0:TT), util(JJ,0:TT), v(-JJ+2:TT)
    real*8 :: l(JJ, 0:TT)
    integer :: tax(0:TT), eps(0:TT)
    logical :: lsra_on, smopec
    integer :: it

    ! initialize variables and government parameters
    call initialize

    ! compute initial long-run equilibrium
    call get_SteadyState

    ! write output
    open(20, file='output.out')
    call output(0, 20)

    ! calculate transition path
    call get_Transition

    do it = 1, TT
        call output(it, 20)
    enddo
    close(20)

    open(21, file='summary.out')
    call output_summary(21)

    ! calculate transition with LSRA payments
    lsra_on = .true.

    if(lsra_on)then
        call get_Transition
        write(21, *)
        call output_summary(21)
    endif
    close(21)

contains


    ! initializes variables and government parameters
    subroutine initialize

        implicit none
        integer :: ij, it, itm

        ! set model parameters
        n_p = 0.2d0
        gy = 0d0
        by = 0d0
        kappa = 0d0
        tauk = 0d0
        eps = 0
        tax = 1
        lsra_on = .false.
        smopec = .false.


        ! human capital profiles for Table 7.1 (uncomment respective line for different rows)
        h = (/1.0d0, 1.0d0, 1.0d0/)
        !h = (/0.5d0, 1.0d0, 1.5d0/)
        !h = (/1.5d0, 1.0d0, 0.5d0/)
        !h = (/0.5d0, 2.0d0, 0.5d0/)


        ! setup for Table 7.2
        !h = (/2.0d0, 2.0d0, 0.0d0/)
        !gy = 0.195d0
        !tax(1:TT) = 3


        ! setup for Table 7.3
        !h = (/2.0d0, 2.0d0, 0.0d0/)
        !gy = 0.195d0
        !kappa(1:TT) = 0.5d0


        ! for everything else you might wanna do
        !by(1:TT) = 0d0
        !eps(1:TT) = 0d0
        !n_p(1:TT) = 0.2d0

        ! initialize tax rates shadow wages and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0
        mu = 0d0
        if(nu > 0d0)mu(JR:JJ,:) = 0.5d0

        ! initialize assets, LSRA payments and debt holdings
        a = 0d0
        c = 0d0
        l = 0d0
        v = 0d0
        YY = 0d0
        BA = 0d0
        BF = 0d0
        TB = 0d0
        TXR = 0d0

        ! size of cohorts in specific year
        do it = 0, TT
            m(1, it) = 1d0
            itm = year(it, 2, 1)
            do ij = 2, JJ
                m(ij, it) = m(ij-1, itm)/(1d0 + n_p(it))
            enddo
        enddo

    end subroutine


    ! solves initial steady state using Gauss-Seidel
    subroutine get_SteadyState

        implicit none
        integer :: iter

        ! initial guess for capital
        KK(0) = 1d0
        LL(0) = 1d0

        do iter = 1, itermax

            ! get prices, decisions, quantities and taxes
            call factor_prices(0)
            call decisions(0)
            call quantities(0)
            call government(0)

            ! check for the number of markets in equilibrium
            if(abs(YY(0) - CC(0) - II(0) - GG(0))/YY(0) < tol)exit
        enddo

        if(iter < itermax)then
            write(*,'(a,i4,a,f16.10)')'Iteration: ', iter, ' Diff: ', &
                                     abs(YY(0)-CC(0)-II(0)-GG(0))/YY(0)
            write(*,*)
        else
            write(*, '(/a/)')'!!! No equilibrium found !!!'
        endif

    end subroutine


    ! solves for transition path using Gauss-Seidel
    subroutine get_Transition

        implicit none
        integer :: it, nmarket, iter

        ! initialize values from initial equilibrium
        a(:, 1) = a(:, 0)
        KK(:) = KK(0)
        LL(:) = LL(0)
        nmarket = 0

        do iter = 1, itermax

            ! get prices, decisions and quantities
            do it = 1, TT
                call factor_prices(it)
            enddo

            do it = 1, TT
                call decisions(it)
            enddo

            if(lsra_on)call lsra()

            do it = 1, TT
                call quantities(it)
            enddo

            do it = 1, TT
                call government(it)
            enddo

            ! check for the number of markets in equilibrium
            nmarket = 0
            do it = 1, TT
                if(abs(YY(it) - CC(it) - II(it) - GG(it) - TB(it))/YY(it) < tol*10d0) &
                    nmarket = nmarket + 1
            enddo

            if(nmarket == TT)exit

        enddo

        if(iter > itermax)then
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        else
            write(*,'(a,i4,a,i4,a,f16.10)')'Iteration: ',iter, &
                ' Markets: ', nmarket,' Diff: ', &
                maxval(abs(YY(:) - CC(:) - II(:) - GG(:) - TB(:))/YY(:))
            write(*,*)
        endif

    end subroutine


    ! subroutine for calculating factor prices in a certain year
    subroutine factor_prices(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij

        ! factor prices and pension payments in year t
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = (1d0-tauk(it))/(1d0-eps(it)*tauk(it))*     &
                    (alpha*(KK(it)/LL(it))**(alpha-1d0)-delta)
        endif
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha

        ! derive after tax prices (incl. shadow wage)
        do ij = 1, JJ
            wn(ij, it) = (h(ij)*w(it) + mu(ij, it))*(1d0-tauw(it)-taup(it))
        enddo
        Rn(it) = 1d0 + r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij

        ! solve cohort that just entered the economy
        c(1, it) = get_Psi(1, it)*get_W(1, it)
        call get_path(1, it)
        if(nu > 0d0)call shadw(1, it)

        ! derive behavior for all other cohorts in year 1 of transition
        if(it == 1)then
            do ij = 2, JJ
                c(ij, it) = get_Psi(ij, it)*get_W(ij, it)
                call get_path(ij, it)
                if(nu > 0d0)call shadw(ij, it)
            enddo
        endif

    end subroutine


    ! present value of resources for household aged ij in year it
    function get_W(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: itp, ijp
        real*8 :: get_W, PRn

        ! get current value of resources
        get_W = wn(ij, it) + pen(ij, it)
        if(it == 1 .and. ij > 1)get_W = get_W + Rn(it)*a(ij, it) + v(-ij+2)
        if(it >= 1 .and. ij == 1)get_W = get_W + v(it)

        ! iterate over remainder of life span
        PRn = 1d0
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_W = get_W + (wn(ijp, itp) + pen(ijp, itp))/PRn
        enddo

    end function


    ! marginal consumption for household aged ij in year it
    function get_Psi(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: itp, itp1, ijp
        real*8 :: get_Psi, PRn, vv(JJ)

        get_Psi = 0d0
        PRn = 1d0

        do ijp = ij, JJ
            itp = year(it, ij, ijp)
            itp1 = year(it, ij, ijp+1)
            vv(ijp) = (1d0 + nu**rho*(wn(ijp, itp)/p(itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            get_Psi = get_Psi + beta**((ijp-ij)*gamma)*(p(itp)/PRn/p(it))**(1d0-gamma)* &
                          vv(ijp)**((1d0-gamma)/(rho-gamma))
            PRn = PRn*Rn(itp1)
        enddo
        get_Psi = vv(ij)/p(it)/get_Psi

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: ijp, itp, itm
        real*8 :: PRn, vv(JJ)

        ! get leisure at age ij
        PRn = 1d0
        vv(ij) = (1d0 + nu**rho*(wn(ij, it)/p(it))**(1d0-rho))**((rho-gamma)/(1d0-rho))
        if(nu > 0d0)l(ij, it) = (wn(ij, it)/nu/p(it))**(-rho)*c(ij, it)

        ! determine consumption and leisure path for remainder of lifetime
        do ijp = ij+1, JJ

            ! get future and previous year as well as interest factor
            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp-1)
            PRn = PRn*Rn(itp)

            ! get consumption and leisure
            vv(ijp) = (1d0 + nu**rho*(wn(ijp, itp)/p(itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            c(ijp, itp) = vv(ijp)/vv(ij)*(beta**(ijp-ij)*PRn*p(it)/p(itp))**gamma*c(ij, it)
            if(nu > 0d0)l(ijp, itp) = (wn(ijp, itp)/nu/p(itp))**(-rho)*c(ijp, itp)

            ! calculate assets
            a(ijp, itp) = wn(ijp-1, itm)*(1d0-l(ijp-1, itm)) + pen(ijp-1, itm) +  &
                       Rn(itm)*a(ijp-1, itm) - p(itm)*c(ijp-1, itm)
            if(itp == 2)a(ijp, itp) = a(ijp, itp) + v(-ijp+3)
            if(itp > 2 .and. ijp == 2)a(ijp, itp) = a(ijp, itp) + v(itm)
        enddo

    end subroutine


    ! compute shadow wages
    subroutine shadw(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: ijp, itp
        real*8 :: mu_new
        logical :: do_all = .false.

        ! iterate over all ages
        do ijp = ij, JJ
            itp = year(it, ij, ijp)

            ! check whether leisure satisfies constraint
            if(l(ijp, itp) > 1d0 .or. do_all)then
                mu_new = c(ijp, itp)**(1d0/rho)*nu*p(itp)/ &
                            (1d0-tauw(itp)-taup(itp)) - h(ijp)*w(itp)
                mu(ijp, itp) = max((1d0-damp)*mu(ijp, itp) + damp*mu_new, 0d0)

                ! once leisure doesn't satisfy constraint, adjust all future shadow wages
                do_all = .true.
            else
                mu(ijp, itp) = (1d0-damp)*mu(ijp, itp)
            endif
        enddo

    end subroutine


    ! subroutine for calculating lsra transfers
    subroutine lsra()

        implicit none
        integer :: it, itp, ij
        real*8 :: sum1, sum2, ustar, PVV

        ! calculate utility for each generation
        do it = 1, TT
            call utility(it)
        enddo

        ! transfers to old generations
        BA = 0d0
        do ij = 2, JJ
            v(-ij+2) = v(-ij+2) + get_W(ij, 1)*((util(ij, 0)/util(ij, 1))**(1d0/egam)-1d0)
            BA(2) = BA(2) + v(-ij+2)*m(ij, 1)
        enddo

        ! long run equilibrium
        PVV  = v(TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum1 = get_W(1, TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum2 = get_W(1, TT)*(util(1, TT)*egam)**(-1d0/egam)*(1d0+r(TT))/(r(TT)-n_p(TT))

        ! transition path
        do it = TT-1, 1, -1
            itp = year(it, 1, 2)
            PVV = PVV*(1d0+n_p(itp))/(1d0+r(itp)) + v(it)
            sum1 = sum1*(1d0+n_p(itp))/(1d0+r(itp)) + get_W(1, it)
            sum2 = sum2*(1d0+n_p(itp))/(1d0+r(itp)) + get_W(1, it)*(util(1,it)*egam)**(-1d0/egam)
        enddo

        ! calculate ustar for future generations
        ustar = ((sum1 - BA(2) - PVV)/sum2)**egam/egam

        ! calculate transfers to future generations and debt of LSRA
        do it = 1, TT
            v(it) = v(it) + get_W(1, it)*((ustar/util(1, it))**(1d0/egam)-1d0)
            if(it == 2)BA(2) = (BA(2) + v(1))/(1d0+n_p(2))
            if(it > 2)BA(it) = ((1d0+r(it-1))*BA(it-1) + v(it-1))/(1d0+n_p(it))
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp, ij

        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        if(it == 0)then
            GG(it) = gy*YY(it)
        else
            GG(it) = GG(0)
        endif

        ! aggregate individual decisions
        CC(it) = 0d0
        AA(it) = 0d0
        LL(it) = 0d0
        do ij = 1, JJ
            CC(it) = CC(it) + c(ij, it)*m(ij, it)
            AA(it) = AA(it) + a(ij, it)*m(ij, it)
            LL(it) = LL(it) + h(ij)*(1d0-l(ij, it))*m(ij, it)
        enddo
        YY(it) = KK(it)**alpha*LL(it)**(1d0-alpha)
        BB(it) = by(itm)*YY(it)

        ! derive capital in small open or closed economy
        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)*(1d0-eps(it)*tauk(it))/(1d0-tauk(it)) + &
                         delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BA(it) - BB(it)
            TB(it) = (1d0+n_p(itp))*BF(itp) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)-BA(it)) + (1d0-damp)*KK(it)
        endif
        II(it) = (1d0+n_p(itp))*KK(itp) - (1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: itp, ij
        real*8 :: taxrev(4)

        itp = year(it, 1, 2)

        taxrev(1) = tauc(it)*CC(it)
        taxrev(2) = tauw(it)*w(it)*LL(it)
        taxrev(3) = taur(it)*r(it)*AA(it)
        taxrev(4) = tauk(it)*(YY(it) - w(it)*LL(it) - (delta+eps(it)*r(it))*KK(it))

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = ((1d0+r(it))*BB(it) + GG(it) - (taxrev(2) + taxrev(3) + &
                           taxrev(4) + (1d0+n_p(itp))*BB(itp)))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = ((1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + taxrev(4)+&
                           (1d0+n_p(itp))*BB(itp)))/(w(it)*LL(it) + r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = ((1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + taxrev(3) + &
                           taxrev(4) + (1d0+n_p(itp))*BB(itp)))/(w(it)*LL(it))
        else
            taur(it) = ((1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + taxrev(2) + &
                           taxrev(4) + (1d0+n_p(itp))*BB(itp)))/(r(it)*AA(it))
        endif
        TXR(it) = sum(taxrev)

        ! get budget balancing social security contribution
        pen(JR:JJ, it) = kappa(it)*w(it)
        Tpen(it) = 0d0
        do ij = JR, JJ
            Tpen(it) = Tpen(it) + pen(ij, it)*m(ij, it)
        enddo
        taup(it) = Tpen(it)/w(it)/LL(it)

    end subroutine


    ! calculates year at which age ij agent is ijp
    function year(it, ij, ijp)

        implicit none
        integer, intent(in) :: it, ij, ijp
        integer :: year

        year = it + ijp - ij

        if(it == 0 .or. year <= 0)year = 0
        if(it == TT .or. year >= TT)year = TT

    end function


    ! subroutine for writing output to a file
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff
        integer :: itp, ij

        ! household utility
        call utility(it)

        itp = year(it, 1, 2)

        ! difference on goods market
        diff = YY(it) - CC(it) - GG(it) - II(it) - TB(it)

        write(fname,'(a,i3/)')'Equilibrium: Year', it
        write(fname,'(a)')' Goods Market  '
        write(fname,'(a)')'      Y      C      G      I     TB      DIFF'
        write(fname,'(5f7.2,f10.4)')YY(it), CC(it), GG(it), II(it), TB(it), diff
        write(fname,'(5f7.2/)')(/YY(it), CC(it), GG(it), II(it), TB(it)/)/YY(it)*100d0
        write(fname,'(a)')' Capital Market  '
        write(fname,'(a)')'      A      K     BB     BA     BF       r'
        write(fname,'(5f7.2,f8.2)')AA(it), KK(it), BB(it), BA(it), BF(it), r(it)
        write(fname,'(5f7.2/)')(/AA(it), KK(it), BB(it), BA(it), BF(it)/)/YY(it)*100d0
        write(fname,'(a)')' Labor Market  '
        write(fname,'(a)')'     LL      w   util'
        write(fname,'(3f7.2/)')LL(it), w(it), util(1, it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk    TXR     DD     rB  '
        write(fname,'(8f7.2/)')tauc(it), tauw(it), taur(it), taup(it), tauk(it), &
                        (/TXR(it), (1d0+n_p(itp))*BB(itp)-BB(it), r(it)*BB(it)/)/YY(it)*100d0

        ! write individual output
        write(fname,'(a)')' Age    cons     leis     wn       mu       pen      asset    Diff      '
        do ij = 1, JJ
            itp = year(it, 1, 2)
            if(ij < JJ)then
                diff = Rn(it)*a(ij, it) + wn(ij, it)*(1d0-l(ij, it)) + pen(ij, it) &
                       - a(ij+1, itp) - p(it)*c(ij, it)
            else
                diff = Rn(it)*a(ij, it) + wn(ij, it)*(1d0-l(ij, it)) + pen(ij, it) &
                       - p(it)*c(ij, it)
            endif
            write(fname,'(i3,7f9.2)')ij, c(ij, it), l(ij, it), wn(ij, it),&
                    mu(ij, it), pen(ij, it), a(ij, it), diff
        enddo
        write(fname,'(/)')

    end subroutine


    ! subroutine for writing output summary to a file
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it, ij
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'         C      A      K      L      r', &
            '      w      Y      B     BA   tauc   tauw   taur', &
            '   taup   tauk      v    HEV   y-d'
        diff = YY(0) - CC(0) - GG(0) - II(0)

        ! current generations
        do ij = JJ, 3, -1
            write(fname,'(i3,98x,2f7.2)')-ij+2, v(-ij+2), &
                  ((util(ij,1)/util(ij,0))**(1d0/egam)-1d0)*100d0
        enddo
        write(fname,'(i3,17f7.2)')0, CC(0), AA(0), KK(0), LL(0), r(0), w(0), &
            YY(0), BB(0), BA(0), tauc(0), tauw(0), taur(0), taup(0), tauk(0), v(0), &
            ((util(2,1)/util(2,0))**(1d0/egam)-1d0)*100d0, diff

        ! future generations
        do it = 1, TT
            diff = YY(it) - CC(it) - GG(it) - II(it)
            write(fname,'(i3,17f7.2)')it, CC(it), AA(it), KK(it), LL(it), r(it), w(it), &
                  YY(it), BB(it), BA(it), tauc(it), tauw(it), taur(it), taup(it), tauk(it), v(it), &
                  ((util(1,it)/util(1,0))**(1d0/egam)-1d0)*100d0, diff
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  itp, ij, ijp

        ! for first generation
        util(1, it) = 0d0
        do ij = 1, JJ
            itp = year(it, 1, ij)
            util(1, it) = util(1, it) + beta**(ij-1)* &
                (c(ij, itp)**erho + nu*l(ij, itp)**erho)**(egam/erho)/egam
        enddo

        ! for current total popilation if year = 0 or 1
        if(it < 2)then
            do ij = 2, JJ
                util(ij, it) = 0d0
                do ijp = ij, JJ
                    itp = year(it, ij, ijp)
                    util(ij, it) = util(ij, it) + beta**(ijp-ij)* &
                        (c(ijp, itp)**erho + nu*l(ijp, itp)**erho)**(egam/erho)/egam
                enddo
            enddo
        endif

    end subroutine

end program
