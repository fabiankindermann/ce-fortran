!##############################################################################
! PROGRAM SKILL_OLG
!
! ##  Two skill class economy
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
program SKILL_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 24    ! 160
    integer, parameter :: JJ      = 3       !  55
    integer, parameter :: JR      = 3       !  41
    integer, parameter :: SS      = 2
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  nu      = 1.5d0
    real*8, parameter ::  rho     = 0.6d0
    real*8, parameter ::  erho    = 1d0 - 1d0/rho
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.2d0
    integer, parameter :: itermax = 200

    ! model variables
    real*8 :: w(0:TT), r(0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT)
    real*8 :: by(0:TT), kappa(0:TT), n_p(0:TT), gy, tauk(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BA(0:TT), BF(0:TT), TB(0:TT)
    real*8 :: Tpen(0:TT), TXR(0:TT)
    real*8 :: mu(JJ, SS, 0:TT), wn(JJ, SS, 0:TT),h(JJ, SS), m(JJ, SS, 0:TT)
    real*8 :: a(JJ, SS, 0:TT), c(JJ, SS, 0:TT), util(JJ, SS, 0:TT), v(SS,-JJ+2:TT)
    real*8 :: l(JJ, SS, 0:TT), pen(JJ, SS, 0:TT)
    integer :: tax(0:TT), eps(0:TT)
    logical :: lsra_on, smopec
    integer :: it

    ! initialize variables and government parameters
    call initialize

    ! compute initial long run equilibrium
    call get_SteadyState

    ! write output
    open(20, file='output.out')
    call output(0, 20)

    ! initialize transitional values
    call get_Transition

    do it = 1, TT
        call output(it, 20)
    enddo
    close(20)

    open(21, file='summary.out')
    call output_summary(21)

    ! get lsra run
    lsra_on = .true.

    if(lsra_on)then
        call get_Transition
        write(21,*)
        call output_summary(21)
    endif
    close(21)

contains


    ! initializes variables and government parameters
    subroutine initialize

        implicit none
        integer :: ij, ik, it, itm

        ! set model parameters
        n_p(0) = 0.2d0
        n_p(1:TT) = 0.2d0
        gy = 0.195d0
!       gy = 0d0
        lsra_on = .false.
        smopec = .false.

        ! set reform values
        by(0:TT) = 0d0
        kappa(0) = 0.0d0
        kappa(1:TT) = 0.5d0
        eps(0) = 0
        eps(1:TT) = 0
        tauk(0) = 0d0
        tauk(1:TT) = 0d0
        tax(0) = 1
        tax(1:TT) = 1

        ! initialize tax rates shadow wages and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0
        mu = 0d0
        if(nu > 0)mu(JJ, :, :) = 0.5

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

        ! human capital profile
        h(:, 1) = 1d0
        h(:, 2) = 2d0
        h(JR:JJ, :) = 0d0

        ! size of cohorts in specific year
        do it = 0, TT
            m(1, 1, it) = 0.0001d0
            m(1, 2, it) = 1d0-m(1, 1, it)
            itm = year(it, 2, 1)
            do ij = 2, JJ
                do ik = 1, SS
                    m(ij,ik,it) = m(ij-1,ik,itm)/(1+n_p(it))
                enddo
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
            if(abs(YY(0) - CC(0) - II(0) - GG(0))/YY(0) < tol) exit
        enddo

        if(iter < itermax)then
            write(*,'(a,i4,a,f16.10)')'Iteration: ',iter, ' Diff: ', &
                                  abs(YY(0)-CC(0)-II(0)-GG(0))/YY(0)
            write(*,*)
        else
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        endif

    end subroutine


    ! solves for transition path using Gauss-Seidel
    subroutine get_Transition()

        implicit none
        integer :: it, nmarket, iter

        ! initialize values from initial equilibrium
        a(:,:,1) = a(:,:,0)

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
        integer :: ik, ij

        ! factor prices and pension payments in year t
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = (1d0-tauk(it))/(1d0-eps(it)*tauk(it))*     &
                    (alpha*(KK(it)/LL(it))**(alpha-1d0)-delta)
        endif
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha
        do ij = 1, JJ
            do ik = 1, SS
                wn(ij, ik, it) = (h(ij, ik)*w(it)+mu(ij, ik, it))*(1d0-tauw(it)-taup(it))
            enddo
        enddo
        Rn(it) = 1d0+r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none
        integer, intent(in) :: it
        integer :: ik, ij

        do ik = 1, SS
            c(1, ik, it) = get_Psi(1, ik, it)*get_W(1, ik, it)
            call getpath(1, ik, it)
            if(nu /= 0)call shadw(1, ik, it)
        enddo


        if(it == 1)then
            do ij = 2, JJ
                do ik = 1, SS
                    c(ij, ik, it) = get_Psi(ij, ik, it)*get_W(ij, ik, it)
                    call getpath(ij, ik, it)
                    if(nu /= 0)call shadw(ij, ik, it)
                enddo
            enddo
        endif


    end subroutine

    ! function for calculation present value of resources for household age j in year it
    function get_W(ij, ik, it)

        implicit none
        integer, intent(in) :: it, ij, ik
        integer :: itp, ijp
        real*8 :: get_W, PRn

        ! Compute index for year corresponding to age i
        get_W = wn(ij, ik, it) + pen(ij, ik, it)
        PRn = 1.0d0

        if(it == 1 .and. ij > 1)get_W = get_W+Rn(it)*a(ij,ik,it)+v(ik,-ij+2)
        if(ij == 1)get_W = get_W + v(ik,it)

        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_W = get_W +  (wn(ijp, ik, itp) + pen(ijp, ik, itp))/PRn
        enddo

    end function

    ! function for calculating the marginal consumption in the initial age
    function get_Psi(ij, ik, it)

        implicit none
        integer, intent(in) :: it, ij, ik
        integer :: itp, itm, ijp
        real*8 :: get_Psi, PRn, vv(JJ)

        get_Psi = 0d0
        PRn = 1d0
        do ijp = ij, JJ
            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp+1)
            vv(ijp) = (1d0+nu**rho*(wn(ijp,ik, itp)/p(itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            get_Psi = get_Psi+beta**((ijp-ij)*gamma)*(p(itp)/PRn/p(it))**(1d0-gamma)*&
                      vv(ijp)**((1d0-gamma)/(rho-gamma))
            PRn = PRn*Rn(itm)
        enddo
        get_Psi = vv(ij)/p(it)/get_Psi

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine getpath(ij, ik, it)

        implicit none
        integer, intent(in) :: it, ij, ik
        integer :: ijp, itp, itm
        real*8 :: PRn, vv(JJ)

        PRn = 1d0
        vv(ij) = (1d0+nu**rho*(wn(ij,ik,it)/p(it))**(1d0-rho))**((rho-gamma)/(1d0-rho))
        if(nu /= 0)l(ij, ik, it) = (wn(ij,ik,it)/nu/p(it))**(-rho)*c(ij,ik,it)
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp-1)
            PRn = PRn*Rn(itp)
            vv(ijp) = (1d0+nu**rho*(wn(ijp,ik,itp)/p(itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            c(ijp,ik,itp) = vv(ijp)/vv(ij)*(beta**(ijp-ij)*PRn*p(it)/p(itp))**gamma*c(ij,ik,it)
            if(nu /= 0)l(ijp, ik, itp) = (wn(ijp,ik,itp)/nu/p(itp))**(-rho)*c(ijp,ik,itp)
            a(ijp,ik,itp) = wn(ijp-1,ik,itm)*(1d0-l(ijp-1,ik,itm))+pen(ijp-1,ik,itm)+ &
                        Rn(itm)*a(ijp-1,ik,itm)-p(itm)*c(ijp-1,ik,itm)
            if(itp == 2)a(ijp,ik,itp)=a(ijp,ik,itp)+v(ik,-ijp+3)
            if(itp > 2 .and. ijp == 2)a(ijp,ik,itp) = a(ijp,ik,itp)+v(ik,itm)
        enddo

    end subroutine

    ! compute shadow wages
    subroutine shadw(ij, ik, it)

        implicit none
        integer, intent(in) :: it, ij, ik
        integer :: ijp, ijm, itp, itm

        do ijp = ij, JJ
            itp = year(it, ij, ijp)
            if(l(ijp,ik,itp) > 1d0 - 5d-2)then
                do ijm = ijp, JJ
                    itm = year(itp, ijp, ijm)
                    mu(ijm,ik,itm) = (1-damp)*mu(ijm,ik,itm)+damp*(c(ijm,ik,itm)**(1d0/rho) &
                                     *nu*p(itm)/(1d0-tauw(itm)-taup(itm)) - h(ijm,ik)*w(itm))
                    if(mu(ijm,ik,itm) < 0) mu(ijm,ik,itm) = 0
                enddo
                exit
            endif
            mu(ijp,ik,itp) = (1-damp)*mu(ijp,ik,itp)
        enddo

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


    ! subroutine for calculating lsra transfers
    subroutine lsra()

        implicit none
        integer :: it, itp, ij, ik
        real*8 :: sum1, sum2, ustar, PVV(SS)

        ! calculate utility for each generation
        do it = 1, TT
            call utility(it)
        enddo

        ! transfers to old generations
        BA = 0
        do ij = 2, JJ
            do ik = 1, SS
                v(ik,-ij+2) = v(ik,-ij+2)+get_W(ij,ik,1)*((util(ij,ik,0)/util(ij,ik,1))**(1d0/egam)-1d0)
                BA(2) = BA(2)+v(ik,-ij+2)*m(ij,ik,1)
            enddo
        enddo


        ! long run equilibrium with low-skill class compensated
        PVV = 0d0
        sum1 = 0d0
        sum2 = 0d0
        PVV(SS) = PVV(SS) + m(1,SS,TT)*v(SS,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum1 = sum1 + m(1,SS,TT)*get_W(1, SS,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum2 = sum2 + m(1,SS,TT)*get_W(1,SS,TT)*(util(1,SS,TT)*egam)**(-1d0/egam)*(1d0+r(TT))/ &
                 (r(TT)-n_p(TT))
        v(1,TT) = v(1,TT)+get_W(1,1,TT)*((util(1,1,0)/util(1,1,TT))**(1d0/egam)-1d0)
        PVV(1) = PVV(1) + m(1,1,TT)*v(1,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))

        ! transition path
        do it = TT-1, 1, -1
            itp = year(it, 1, 2)
            PVV(SS) = PVV(SS)*(1d0+n_p(itp))/(1d0+r(itp)) + m(1,SS,it)*v(SS,it)
            sum1 = sum1*(1d0+n_p(itp))/(1d0+r(itp)) + m(1,SS,it)*get_W(1,SS,it)
            sum2 = sum2*(1d0+n_p(itp))/(1d0+r(itp)) + m(1,SS,it)*get_W(1,SS,it)*&
                   (util(1,SS,it)*egam)**(-1d0/egam)
            v(1,it) = v(1,it)+get_W(1,1,it)*((util(1,1,0)/util(1,1,it))**(1d0/egam)-1d0)
            PVV(1) = PVV(1)*(1d0+n_p(itp))/(1d0+r(itp)) + m(1,1,it)*v(1,it)
        enddo

        ! calculate ustar for future generations
        ustar = ((sum1-BA(2)-sum(PVV))/sum2)**egam/egam

        ! calculate transfers to future generations and debt of LSRA
        do it = 1, TT
            if(it > 2)BA(it) = (1d0+r(it-1))*BA(it-1)
            v(SS,it) = v(SS,it) + get_W(1,SS,it)*((ustar/util(1,SS,it))**(1d0/egam)-1d0)
            do ik = 1, SS
                if(it > 1)BA(it) = BA(it) + m(1,ik,it-1)*v(ik,it-1)
            enddo
            BA(it) = BA(it)/(1d0+n_p(it))
        enddo



    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp, ij, ik

        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        if(it == 0)then
            GG(it) = gy*YY(it)
        else
            GG(it) = GG(0)
        endif

        CC(it) = 0d0
        AA(it) = 0d0
        LL(it) = 0d0
        do ij = 1, JJ
            do ik = 1, SS
                CC(it) = CC(it) + c(ij, ik, it)*m(ij, ik, it)
                AA(it) = AA(it) + a(ij, ik, it)*m(ij, ik, it)
                LL(it) = LL(it) + h(ij,ik)*(1d0-l(ij, ik, it))*m(ij, ik, it)
            enddo
        enddo
        YY(it) = KK(it)**alpha * LL(it)**(1d0-alpha)
        BB(it) = by(itm)*YY(it)

        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)*(1d0-eps(it)*tauk(it))/(1d0-tauk(it))+ &
                   delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BA(it) - BB(it)
            TB(it) = (1d0+n_p(itp))*BF(itp) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)- BA(it)) + (1d0-damp)*KK(it)
        endif
        II(it) = (1d0+n_p(itp))*KK(itp)-(1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: itp, ij, ik
        real*8 :: taxrev(4)

        itp = year(it, 1, 2)

        taxrev(1) = tauc(it)*CC(it)
        taxrev(2) = tauw(it)*w(it)*LL(it)
        taxrev(3) = taur(it)*r(it)*AA(it)
        taxrev(4) = tauk(it)*(YY(it)-w(it)*LL(it)-(delta+eps(it)*r(it))*KK(it))

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(2) + taxrev(3) + &
                       taxrev(4) + (1d0+n_p(itp))*BB(itp)))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(4)+&
                       (1d0+n_p(itp))*BB(itp)))/(w(it)*LL(it)+r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(3) + &
                      taxrev(4) + (1d0+n_p(itp))*BB(itp)))/(w(it)*LL(it))
        else
            taur(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(2) + &
                      taxrev(4) + (1d0+n_p(itp))*BB(itp)))/(r(it)*AA(it))
        endif
        TXR(it) = sum(taxrev)


        pen(JR:JJ,:,it) = kappa(it)*w(it)
        ! get budget balancing social security contribution
        Tpen(it) = 0d0
        do ij = JR, JJ
            do ik = 1, SS
                Tpen(it) = Tpen(it) + pen(ij,ik,it)*m(ij,ik,it)
            enddo
        enddo
        taup(it) = Tpen(it)/w(it)/LL(it)

    end subroutine


    ! subroutine for writing output to a file or the console (fname = 0)
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff
        integer :: itp, ij, ik


        ! household utility
        call utility(it)

        itp = year(it, 1, 2)

        ! difference on goods market
        diff = YY(it)-CC(it)-GG(it)-II(it)-TB(it)

        write(fname,'(a,i3/)')'Equilibrium: Year', it
        write(fname,'(a)')' Goods Market  '
        write(fname,'(a)')'      Y      C      G      I     TB      DIFF'
        write(fname,'(5f7.2,f10.4)')YY(it),CC(it),GG(it),II(it),TB(it),diff
        write(fname,'(5f7.2/)')YY(it)/YY(it),CC(it)/YY(it),GG(it)/YY(it),II(it)/YY(it),TB(it)/YY(it)
        write(fname,'(a)')' Capital Market  '
        write(fname,'(a)')'      A      K     BB     BA     BF       r'
        write(fname,'(5f7.2,f8.2)')AA(it), KK(it), BB(it), BA(it), BF(it), r(it)
        write(fname,'(5f7.2/)')AA(it)/YY(it),KK(it)/YY(it),BB(it)/YY(it),BA(it)/YY(it),BF(it)/YY(it)
        write(fname,'(a)')' Labor Market  '
        write(fname,'(a)')'     LL      w   util'
        write(fname,'(4f7.2/)')LL(it),w(it),util(1,1,it), util(1,2,it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk    TXR     DD     rB  '
        write(fname,'(8f7.2/)')tauc(it),tauw(it),taur(it),taup(it),tauk(it),TXR(it)/YY(it),&
                            ((1+n_p(itp))*BB(itp)-BB(it))/YY(it),r(it)*BB(it)/YY(it)

        if(fname > 1)then
            write(fname,'(a)')' Age    cons     leis     wn       mu       pen      asset    Diff      '
            do ik = 1, SS
                do ij = 1, JJ
                    itp = year(it, 1, 2)
                    if(ij < JJ)then
                        diff = Rn(it)*a(ij,ik,it)+wn(ij,ik,it)*(1d0-l(ij,ik,it)) + &
                                  pen(ij,ik,it)-a(ij+1,ik,itp)-p(it)*c(ij,ik,it)
                    else
                        diff = Rn(it)*a(ij,ik,it)+wn(ij,ik,it)*(1d0-l(ij,ik,it)) + &
                                  pen(ij,ik,it)-p(it)*c(ij,ik,it)
                    endif
                    write(fname,'(i3,7f9.2)')ij,c(ij,ik,it),l(ij,ik,it),wn(ij,ik,it), &
                          mu(ij,ik,it),pen(ij,ik,it),a(ij,ik,it),diff
                enddo
            enddo
            write(fname,'(a/)')'                                           '
        endif


    end subroutine


    ! subroutine for writing output summary to a file or the console (fname = 0)
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it, ij, ik
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'         C      A      K      L      r', &
            '      w      Y      B     BA   tauc   tauw   taur', &
            '   taup   tauk      v      v    HEV   HEV     y-d'
        diff = YY(0)-CC(0)-GG(0)-II(0)
        ! current generations
        do ij = JJ, 3, -1
            write(fname,'(i3,98x,4f7.2)')-ij+2, (v(ik,-ij+2),ik = 1, SS), &
                  (((util(ij,ik,1)/util(ij,ik,0))**(1d0/egam)-1d0)*100d0, ik = 1, SS)
            write(fname,'(i3,19f7.2)')0,CC(0),AA(0),KK(0),LL(0),r(0),w(0), &
                  YY(0),BB(0),0.00,tauc(0),tauw(0),taur(0),taup(0),tauk(0),(v(ik,0),ik = 1, SS), &
                 (((util(2,ik,1)/util(2,ik,0))**(1d0/egam)-1d0)*100d0,ik = 1,SS),diff
        enddo

        ! future generations
        do it = 1, TT
            diff = YY(it)-CC(it)-GG(it)-II(it)
            write(fname,'(i3,19f7.2)')it,CC(it),AA(it),KK(it),LL(it),r(it),w(it), &
                YY(it),BB(it),BA(it),tauc(it),tauw(it),taur(it),taup(it),tauk(it),(v(ik,it),ik = 1, SS), &
                (((util(1,ik,it)/util(1,ik,0))**(1d0/egam)-1d0)*100d0,ik = 1, SS),diff
        enddo


    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  ik, ij, ijp, itp

        util(1,:,it) = 0d0
        do ij = 1, JJ
            do ik = 1, SS
                itp = year(it, 1, ij)
                util(1,ik,it) = util(1,ik,it)+beta**(ij-1)*(c(ij,ik,itp)**erho+nu*l(ij,ik,itp)**erho)**(egam/erho)/egam
            enddo
        enddo

        if(it < 2)then
            do ij = 2, JJ
                do ik = 1, SS
                    util(ij,ik,it) = 0d0
                    do ijp = ij, JJ
                        itp = year(it, ij, ijp)
                        util(ij,ik,it) = util(ij,ik,it)+beta**(ijp-ij)*(c(ijp,ik,itp)**erho+nu*l(ijp,ik,itp)**erho) &
                                                                       **(egam/erho)/egam
                    enddo
                enddo
            enddo
        endif


    end subroutine

end program
