!##############################################################################
! PROGRAM TWO_COUNTRY_OLG
!
! ##  Two country OLG model
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
program TWO_COUNTRY_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 24    ! 160
    integer, parameter :: JJ      = 3       !  55
    integer, parameter :: JR      = 3       ! 41
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
    real*8, parameter ::  damp    = 0.1d0
    integer, parameter :: itermax = 200

    ! model variables
    real*8 :: w(SS,0:TT), r(0:TT), mu(SS,JJ, 0:TT), wn(SS,JJ, 0:TT), Rn(SS,0:TT), p(SS,0:TT)
    real*8 :: tauw(SS,0:TT), taur(SS,0:TT), tauc(SS,0:TT), taup(SS,0:TT), pen(SS,JJ,0:TT)
    real*8 :: by(SS,0:TT), kappa(SS,0:TT), n_p(0:TT), gy(SS), tauk(SS,0:TT)
    real*8 :: KK(SS,0:TT), LL(SS,0:TT), YY(SS,0:TT), AA(SS,0:TT), CC(SS,0:TT), II(SS,0:TT)
    real*8 :: BB(SS,0:TT), GG(SS,0:TT), BA(SS,0:TT), TB(SS,0:TT)
    real*8 :: Tpen(SS,0:TT), TXR(SS,0:TT), h(SS,JJ), M(SS,JJ,0:TT)
    real*8 :: a(SS,JJ,0:TT), c(SS,JJ,0:TT), util(SS,JJ,0:TT), v(SS,-JJ+2:TT)
    real*8 :: l(SS,JJ, 0:TT)
    integer :: tax(SS,0:TT), eps(SS,0:TT)
    logical :: lsra_on
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
        integer :: ij, k, it, itm

        ! set model parameters
        n_p(0:TT) = 0.2d0
        gy(:) = 0.195d0

        lsra_on = .false.

        ! set reform values
        by(:,0:TT) = 0d0
        kappa(:,0) = 0.0d0
        kappa(1,1:TT) = 0d0
        kappa(2,1:TT) = 0.0d0
        eps(:,0) = 0
        eps(:,1:TT) = 0
        tauk(:,0) = 0d0
        tauk(:,1:TT) = 0d0
        tax(1,:) = 1
        tax(2,0) = 1
        tax(2,1:TT) = 2

        ! initialize tax rates shadow wages and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0
        mu = 0d0
        if(nu > 0) mu(:,JJ,:) = 0.5

        ! initialize assets, LSRA payments and debt holdings
        a = 0d0
        c = 0d0
        l = 0d0
        v = 0d0
        YY = 0d0
        BA = 0d0
        TB = 0d0
        TXR = 0d0

        ! human capital profile
        do ij = 1, JR-1
!             h(ij) = 1.0d0 + 0.05*ij-0.001*ij**2
            h(:,ij) = 2.0d0
        enddo
        h(:,JR:JJ) = 0d0

        ! size of cohorts in specific year
        do it = 0, TT
            m(1,1,it) = 0.5d0
            m(2,1,it) = 1d0 - m(1,1,it)
            itm = year(it, 2, 1)
            do k = 1, 2
                do ij = 2, JJ
                    m(k,ij,it) = m(k,ij-1,itm)/(1+n_p(it))
                enddo
            enddo
        enddo

    end subroutine


    ! solves initial steady state using Gauss-Seidel
    subroutine get_SteadyState

        implicit none
        integer :: iter

        ! initial guess for capital
        KK(:,0) = 1d0
        LL(:,0) = 1d0

        do iter = 1, itermax

            ! get prices, decisions, quantities and taxes
            call factor_prices(0)
            call decisions(0)
            call quantities(0)
            call government(0)

            ! check for the number of markets in equilibrium
            if(abs(YY(1,0) - CC(1,0) - II(1,0) - GG(1,0)-TB(1,0))/YY(1,0) < tol .and. abs(YY(2,0) &
                     - CC(2,0) - II(2,0) - GG(2,0)-TB(2,0))/YY(2,0) < tol)exit
        enddo

        if(iter < itermax)then
            write(*,'(a,i4,a,f10.7)')'Iteration: ',iter, &
                                    ' Country 1: ', abs(YY(1,0)-CC(1,0)-II(1,0)-GG(1,0)-TB(1,0))/YY(1,0)
            write(*,'(15x,a,f10.7)')' Country 2: ', abs(YY(2,0)-CC(2,0)-II(2,0)-GG(2,0)-TB(2,0))/YY(2,0)
            write(*,*)
        else
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        endif

    end subroutine


    ! solves for transition path using Gauss-Seidel
    subroutine get_Transition()

        implicit none
        integer :: it, nmarket, iter!, iterh

        ! initialize values from initial equilibrium
        a(:,:,1) = a(:,:,0)

        KK(1,:) = KK(1,0)
        KK(2,:) = KK(2,0)
        LL(1,:) = LL(1,0)
        LL(2,:) = LL(2,0)
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
                if(abs(YY(1,it) - CC(1,it) - II(1,it) - GG(1,it) - TB(1,it))/YY(1,it) < tol*10d0 .and. abs(YY(2,it)&
                        - CC(2,it) - II(2,it) - GG(2,it) - TB(2,it))/YY(2,it) < tol*10d0)nmarket = nmarket + 1
            enddo

            if(nmarket == TT)exit

        enddo

        if(iter > itermax)then
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        else
            write(*,'(a,i4,a,i4,a,f12.6)')'Iteration: ',iter, ' Markets: ', nmarket,&
                    ' Country 1: ', maxval(abs(YY(1,:) - CC(1,:) - II(1,:) - GG(1,:) - TB(1,:))/YY(1,:))
            write(*,'(29x,a,f12.6)')' Country 2: ', maxval(abs(YY(2,:) - CC(2,:) - II(2,:) - GG(2,:) - TB(2,:))/YY(2,:))
            write(*,*)
        endif

    end subroutine


    ! subroutine for calculating factor prices in a certain year
    subroutine factor_prices(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, ik

        ! factor and consumer prices in year t
        r(it) = (1d0-tauk(1,it))/(1d0-eps(1,it)*tauk(1,it))*     &
                    (alpha*(KK(1,it)/LL(1,it))**(alpha-1d0)-delta)

        do ik = 1, SS
            w(ik,it) = (1d0-alpha)*(KK(ik,it)/LL(ik,it))**alpha
            do ij = 1, JJ
                wn(ik, ij, it) = (h(ik,ij)*w(ik,it)+mu(ik,ij,it))*(1d0-tauw(ik,it)-taup(ik,it))
            enddo
            Rn(ik,it) = 1d0+r(it)*(1d0-taur(ik,it))
            p(ik,it) = 1d0 + tauc(ik,it)
        enddo

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it
        integer :: ik, ij

        do ik = 1, SS
            c(ik, 1, it) = get_Psi(ik, 1, it)*get_W(ik, 1, it)
            call get_path(ik, 1, it)
            if(nu /= 0)call shadw(ik, 1, it)
        enddo


        if(it == 1)then
            do ik = 1, SS
                do ij = 2, JJ
                    c(ik, ij, it) = get_Psi(ik,ij,it)*get_W(ik,ij,it)
                    call get_path(ik, ij, it)
                    if(nu /= 0)call shadw(ik,ij,it)
                enddo
            enddo
        endif

    end subroutine


    ! function for calculation present value of resources for household age ijp in year ij
    function get_W(ik, ij, it)

        implicit none

        ! variable declaration
        integer, intent(in) :: ij, it, ik
        integer :: itp, ijp
        real*8 :: get_W, PRn

        ! Compute index for year corresponding to age it
        get_W = wn(ik, ij, it) + pen(ik, ij, it)
        PRn = 1.0d0

        if(it == 1 .and. ij > 1) get_W = get_W+Rn(ik,it)*a(ik,ij,it)+v(ik,-ij+2)
        if(ij == 1)  get_W = get_W+v(ik,it)

        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(ik,itp)
            get_W = get_W +  (wn(ik, ijp, itp) + pen(ik, ijp, itp))/PRn
        enddo

    end function


    ! function for calculating the marginal consumption in the initial age
    function get_Psi(ik, ij, it)

        implicit none

        ! variable declaration
        integer, intent(in) :: ij, it, ik
        integer :: itp, itp1, ijp
        real*8 :: get_Psi, PRn, vv(JJ)

        get_Psi = 0d0
        PRn = 1d0
        do ijp = ij, JJ
            itp = year(it, ij, ijp)
            itp1 = year(it, ij, ijp+1)
            vv(ijp) = (1d0+nu**rho*(wn(ik,ijp,itp)/p(ik,itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            get_Psi = get_Psi+beta**((ijp-ij)*gamma)*(p(ik,itp)/PRn/p(ik,it))**(1d0-gamma)*&
                        vv(ijp)**((1d0-gamma)/(rho-gamma))
            PRn = PRn*Rn(ik,itp1)
        enddo
        get_Psi = vv(ij)/p(ik,it)/get_Psi

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ik, ij, it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij, ik
        integer :: ijp, itp, itm
        real*8 :: PRn, vv(JJ)

        PRn = 1d0
        vv(ij) = (1d0+nu**rho*(wn(ik,ij,it)/p(ik,it))**(1d0-rho))**((rho-gamma)/(1d0-rho))
        if(nu /= 0)l(ik, ij, it) = (wn(ik,ij,it)/nu/p(ik,it))**(-rho)*c(ik,ij,it)
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp-1)
            PRn = PRn*Rn(ik,itp)
            vv(ijp) = (1d0+nu**rho*(wn(ik,ijp,itp)/p(ik,itp))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            c(ik, ijp, itp) = vv(ijp)/vv(ij)*(beta**(ijp-ij)*PRn*p(ik,it)/p(ik,itp))**gamma*c(ik,ij,it)
            if(nu /= 0)l(ik, ijp, itp) = (wn(ik,ijp,itp)/nu/p(ik,itp))**(-rho)*c(ik,ijp,itp)
            a(ik, ijp, itp) = wn(ik,ijp-1,itm)*(1d0-l(ik,ijp-1,itm))+pen(ik,ijp-1,itm)+ &
                            Rn(ik,itm)*a(ik,ijp-1,itm)-p(ik,itm)*c(ik,ijp-1,itm)
            if(itp == 2)a(ik,ijp,itp) = a(ik, ijp,itp)+v(ik,-ijp+3)
            if(itp > 2 .and. ijp == 2)a(ik,ijp,itp) = a(ik,ijp,itp)+v(ik,itm)
        enddo

    end subroutine


    ! compute shadow wages
    subroutine shadw(ik, ij, it)

        implicit none

        integer, intent(in) :: it, ij, ik
        integer :: ijp, ijp2, itp, itp2

        do ijp = ij, JJ
            itp = year(it, ij, ijp)
            if(l(ik, ijp, itp) > 1d0 - 5d-2)then
                do ijp2 = ijp, JJ
                    itp2 = year(itp, ijp, ijp2)
                    mu(ik, ijp2, itp2) = (1-damp)*mu(ik,ijp2, itp2)+damp*(c(ik,ijp2,itp2)**(1d0/rho) &
                                *nu*p(ik,itp2)/(1d0-tauw(ik,itp2)-taup(ik,itp2)) - h(ik,ijp2)*w(ik,itp2))
                    if(mu(ik, ijp2, itp2) < 0) mu(ik, ijp2, itp2) = 0
                enddo
                exit
            endif
            mu(ik, ijp, itp) = (1-damp)*mu(ik, ijp, itp)
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
                v(ik,-ij+2) = v(ik,-ij+2)+get_W(ik,ij,1)*((util(ik,ij,0)/util(ik,ij,1))**(1d0/egam)-1d0)
                BA(ik,2) = BA(ik,2)+v(ik,-ij+2)*m(ik,ij,1)
            enddo
        enddo

        ! long run equilibrium
        PVV = 0d0
        sum1 = 0d0
        sum2 = 0d0
        PVV(SS) = PVV(SS) + m(SS,1,TT)*v(SS,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum1 = sum1 + m(SS,1,TT)*get_W(SS,1,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum2 = sum2 + m(SS,1,TT)*get_W(SS,1,TT)*(util(SS,1,TT)*egam)**(-1d0/egam)*(1d0+r(TT))/ &
                    (r(TT)-n_p(TT))
        v(1,TT) = v(1,TT)+get_W(1,1,TT)*((util(1,1,0)/util(1,1,TT))**(1d0/egam)-1d0)
        PVV(1) = PVV(1) + m(1,1,TT)*v(1,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))

        ! transition path
        do it = TT-1, 1, -1
            itp = year(it, 1, 2)
            PVV(SS) = PVV(SS)*(1d0+n_p(itp))/(1d0+r(itp)) + m(SS,1,it)*v(SS,it)
            sum1 = sum1*(1d0+n_p(itp))/(1d0+r(itp)) + m(SS,1,it)*get_W(SS,1,it)
            sum2 = sum2*(1d0+n_p(itp))/(1d0+r(itp)) + m(SS,1,it)*get_W(SS,1,it)*&
                    (util(SS,1,it)*egam)**(-1d0/egam)
            v(1,it) = v(1,it)+get_W(1,1,it)*((util(1,1,0)/util(1,1,it))**(1d0/egam)-1d0)
            PVV(1) = PVV(1)*(1d0+n_p(itp))/(1d0+r(itp)) + m(1,1,it)*v(1,it)
        enddo

        ! calculate ustar for future generations
        ustar = ((sum1-sum(BA(:,2))-sum(PVV))/sum2)**egam/egam

        ! calculate transfers to future generations and debt of LSRA
        do it = 1, TT
            v(SS,it) = v(SS,it) + get_W(SS,1,it)*((ustar/util(SS,1,it))**(1d0/egam)-1d0)
            do ik = 1, SS
                if(it == 2) BA(ik,2) = (BA(ik,2) + m(ik,1,1)*v(ik,1))/(1d0+n_p(2))
                if(it > 2) BA(ik,it) = ((1d0+r(it-1))*BA(ik,it-1) + m(ik,1,it-1)*v(ik,it-1))/(1d0+n_p(it))
            enddo
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp, ij, ik
        real :: n_ph

        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        if(it == 0)then
            GG(:,it) = gy(:)*YY(:,it)
        else
            GG(:,it) = GG(:,0)
        endif

        CC(:,it) = 0d0
        AA(:,it) = 0d0
        LL(:,it) = 0d0
        do ik = 1, 2
            do ij = 1, JJ
                CC(ik,it) = CC(ik,it)+c(ik,ij,it)*m(ik,ij,it)
                AA(ik,it) = AA(ik,it)+a(ik,ij,it)*m(ik,ij,it)
                LL(ik,it) = LL(ik,it)+h(ik,ij)*(1d0-l(ik,ij, it))*m(ik,ij,it)
            enddo
            YY(ik,it) = KK(ik,it)**alpha*LL(ik,it)**(1d0-alpha)
            BB(ik,it) = by(ik,it)*YY(ik,it)
        enddo

        KK(2,it) = LL(2,it)*((r(it)*(1d0-eps(2,it)*tauk(2,it))/(1d0-tauk(2,it))+ &
                    delta)/alpha)**(1d0/(alpha-1d0))
        KK(1,it) = damp*(sum(AA(:,it))-sum(BB(:,it))-sum(BA(:,it))-KK(2,it))+(1d0-damp)*KK(1,it)

        do ik = 1, SS
            II(ik,it) = (1d0+n_p(itp))*KK(ik,itp)-(1d0-delta)*KK(ik,it)
            if(ik == 1)TB(ik,it) = YY(ik,it) - CC(ik,it)  - II(ik,it) - GG(ik,it)
            if(ik == 2)TB(2,it) = - TB(1,it)
        enddo

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: itp, ij, ik
        real*8 :: taxrev(SS,4)

        itp = year(it, 1, 2)
        do ik = 1, SS
                taxrev(ik,1) = tauc(ik,it)*CC(ik,it)
                taxrev(ik,2) = tauw(ik,it)*w(ik,it)*LL(ik,it)
                taxrev(ik,3) = taur(ik,it)*r(it)*AA(ik,it)
                taxrev(ik,4) = tauk(ik,it)*(YY(ik,it)-w(ik,it)*LL(ik,it)-(delta+eps(ik,it)*r(it))*KK(ik,it))

                ! get budget balancing tax rate
                if(tax(ik,it) == 1)then
                    tauc(ik,it) = ((1d0+r(it))*BB(ik,it)+GG(ik,it)-(taxrev(ik,2) + taxrev(ik,3) + &
                        taxrev(ik,4) + (1d0+n_p(itp))*BB(ik,itp)))/CC(ik,it)
                elseif(tax(ik,it) == 2)then
                    tauw(ik,it) = ((1d0+r(it))*BB(ik,it)+GG(ik,it)-(taxrev(ik,1) + taxrev(ik,4)+&
                        (1d0+n_p(itp))*BB(ik,itp)))/(w(ik,it)*LL(ik,it)+r(it)*AA(ik,it))
                    taur(ik,it) = tauw(ik,it)
                elseif(tax(ik,it) == 3)then
                    tauw(ik,it) = ((1d0+r(it))*BB(ik,it)+GG(ik,it)-(taxrev(ik,1) + taxrev(ik,3) + &
                        taxrev(ik,4) + (1d0+n_p(itp))*BB(ik,itp)))/(w(ik,it)*LL(ik,it))
                else
                    taur(ik,it) = ((1d0+r(it))*BB(ik,it)+GG(ik,it)-(taxrev(ik,1) + taxrev(ik,2) + &
                        taxrev(ik,4) + (1d0+n_p(itp))*BB(ik,itp)))/(r(it)*AA(ik,it))
                endif
                TXR(ik,it) = sum(taxrev(ik,:))


                pen(ik,JR:JJ,it) = kappa(ik,it)*w(ik,it)
                ! get budget balancing social security contribution
                Tpen(ik,it) = 0d0
                do ij = JR, JJ
                    Tpen(ik,it) = Tpen(ik,it) + pen(ik,ij,it)*m(ik,ij,it)
                enddo
                taup(ik,it) = Tpen(ik,it)/w(ik,it)/LL(ik,it)
            enddo

    end subroutine


    ! subroutine for writing output to a file or the console (fname = 0)
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff(2)
        integer :: itp, ij, ik

        ! household utility
        call utility(it)

        itp = year(it, 1, 2)

        do ik = 1, 2
            ! difference on goods market
            diff(ik) = YY(ik,it)-CC(ik,it)-GG(ik,it)-II(ik,it)-TB(ik,it)

            write(fname,'(a,i3, a, i3/)')'Equilibrium: Year', it, ' Land ', ik
            write(fname,'(a)')' Goods Market  '
            write(fname,'(a)')'      Y      C      G      I     TB      DIFF'
            write(fname,'(5f7.4,f10.4)')YY(ik,it),CC(ik,it),GG(ik,it),II(ik,it),TB(ik,it),diff(ik)
            write(fname,'(5f7.2/)')YY(ik,it)/YY(ik,it),CC(ik,it)/YY(ik,it),GG(ik,it)/YY(ik,it),&
                                II(ik,it)/YY(ik,it),TB(ik,it)/YY(ik,it)
            write(fname,'(a)')' Capital Market  '
            write(fname,'(a)')'      A      k     BB     r'
            write(fname,'(3f7.2,f8.2)')AA(ik,it), KK(ik,it), BB(ik,it), r(it)
            write(fname,'(3f7.2/)')AA(ik,it)/YY(ik,it),KK(ik,it)/YY(ik,it),BB(ik,it)/YY(ik,it)
            write(fname,'(a)')' Labor Market  '
            write(fname,'(a)')'     LL      w   util'
            write(fname,'(3f7.2/)')LL(ik,it),w(ik,it),util(ik,1,it)
            write(fname,'(a)')' GOVERMENT  '
            write(fname,'(a)')'   tauc   tauw   taur   taup   tauk    TXR     DD     rB  '
            write(fname,'(8f7.2)')tauc(ik,it),tauw(ik,it),taur(ik,it),taup(ik,it),tauk(ik,it),TXR(ik,it),&
                            ((1+n_p(itp))*BB(ik,itp)-BB(ik,it)),r(it)*BB(ik,it)
            write(fname,'(35x,3f7.2/)')TXR(ik,it)/YY(ik,it),((1+n_p(itp))*BB(ik,itp)-BB(ik,it))/YY(ik,it),&
                            r(it)*BB(ik,it)/YY(ik,it)


            if(fname > 1)then
                write(fname,'(a)')' Age    cons     leis     wn       mu       pen      asset    Diff      '
                do ij = 1, JJ
                    itp = year(it, 1, 2)
                    if(ij < JJ)then
                        diff(ik) = Rn(ik,it)*a(ik,ij,it)+wn(ik,ij,it)*(1d0-l(ik,ij,it))+pen(ik,ij,it)-&
                                    a(ik,ij+1,itp)-p(ik,it)*c(ik,ij,it)
                    else
                        diff(ik) = Rn(ik,it)*a(ik,ij,it)+wn(ik,ij,it)*(1d0-l(ik,ij,it))+pen(ik,ij,it)-&
                                    p(ik,it)*c(ik,ij,it)
                    endif
                    write(fname,'(i3,7f9.2)')ij,c(ik,ij,it),l(ik,ij,it),wn(ik,ij,it),mu(ik,ij,it),&
                                            pen(ik,ij,it),a(ik,ij,it),diff(ik)
                enddo
                write(fname,'(a/)')'                                           '
            endif
        enddo

    end subroutine


    ! subroutine for writing output summary to a file or the console (fname = 0)
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it, ij, ik
        real*8 :: diff(SS)

        do ik = 1, SS
            ! headline
            write(fname,'(/a,a,a)')'         C      A      k      L      r', &
                '      w      Y      B     BA   TB     tauc   tauw   taur', &
                '   taup   tauk      v    HEV   y-d'
            write(fname, '(a, i3)')'Land', ik
                diff(ik) = YY(ik,0)-CC(ik,0)-GG(ik,0)-II(ik,0)-TB(ik,0)
            ! current generations
            do ij = JJ, 3, -1
                write(fname,'(i3,105x,2f7.2)')-ij+2, v(ik,-ij+2), &
                        ((util(ik,ij,1)/util(ik,ij,0))**(1d0/egam)-1d0)*100d0
            enddo
            write(fname,'(i3,17f7.2,f8.5)')0,CC(ik,0),AA(ik,0),KK(ik,0),LL(ik,0),r(0),w(ik,0), &
                    YY(ik,0),BB(ik,0),BA(ik,0),TB(ik,0),tauc(ik,0),tauw(ik,0),taur(ik,0),taup(ik,0), &
                    tauk(ik,0),v(ik,0),((util(ik,2,1)/util(ik,2,0))**(1d0/egam)-1d0)*100d0,diff(ik)

            ! future generations
            do it = 1, TT
                diff(ik) = YY(ik,it)-CC(ik,it)-GG(ik,it)-II(ik,it)-TB(ik,it)
                write(fname,'(i3,17f7.2,f8.5)')it,CC(ik,it),AA(ik,it),KK(ik,it),LL(ik,it),r(it),w(ik,it), &
                        YY(ik,it),BB(ik,it),sum(BA(:,it)),TB(ik,it),tauc(ik,it),tauw(ik,it),taur(ik,it),taup(ik,it), &
                        tauk(ik,it),v(ik,it),((util(ik,1,it)/util(ik,1,0))**(1d0/egam)-1d0)*100d0,diff(ik)
            enddo
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  itp, ij, ijp, ik

        util(:,1,it) = 0d0

        do ik = 1, SS
            do ijp = 1, JJ
                itp = year(it, 1, ijp)
                util(ik,1,it) = util(ik,1,it)+beta**(ijp-1)*(c(ik,ijp,itp)**erho+ &
                                nu*l(ik,ijp,itp)**erho)**(egam/erho)/egam
            enddo
        enddo

        if(it < 2)then
            do ik = 1, SS
                do ij = 2, JJ
                    util(ik,ij,it) = 0d0
                    do ijp = ij, JJ
                    itp = year(it, ij, ijp)
                    util(ik,ij,it) = util(ik,ij,it) + beta**(ijp-ij)*(c(ik,ijp,itp)**erho+nu*l(ik,ijp,itp)**erho) &
                                                                        **(egam/erho)/egam
                    enddo
                enddo
            enddo
        endif

    end subroutine

end program
