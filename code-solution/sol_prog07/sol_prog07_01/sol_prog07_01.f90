!##############################################################################
! PROGRAM TRVL_OLG
!
! ##  CD and CES utility model
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program TRVL_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 24      ! 160
    integer, parameter :: JJ      = 3       !  55
    integer, parameter :: JR      = 3       ! 41
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  nu      = 0.4d0    ! CES: nu=1.5   CD: nu=0.4
    real*8, parameter ::  rho     = 0.6d0
    real*8, parameter ::  erho    = 1d0 - 1d0/rho
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.2d0
    integer, parameter :: itermax = 200
    logical, parameter :: CES     = .false.
    logical, parameter :: CD      = .true.

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
        n_p(0) = 0.2d0
        n_p(1:TT) = 0.2d0
        gy = 0.195d0
!       gy = 0d0
        lsra_on = .false.
        smopec = .false.

        ! set reform values
        by(0:TT) = 0d0
        kappa(0) = 0.0d0
        kappa(1:TT) = 0.0d0
        eps(0) = 0
        eps(1:TT) = 0
        tauk(0) = 0d0
        tauk(1:TT) = 0d0
        tax(0) = 1
        tax(1:TT) = 3

        ! initialize tax rates shadow wages and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0
        mu = 0d0
        if(CES .or. CD) mu(JJ,:) = 0.5

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
        do ij = 1, JR-1
!          h(ij) = 1.0d0 + 0.05*ij-0.001*ij**2
            h(ij) = 2.0d0
        enddo
        h(JR:JJ) = 0d0

        ! size of cohorts in specific year
        do it = 0, TT
            m(1,it) = 1d0
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
    subroutine get_Transition()

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
        do ij = 1, JJ
            wn(ij, it) = (h(ij)*w(it) + mu(ij, it))*(1d0-tauw(it)-taup(it))
        enddo
        Rn(it) = 1d0 + r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it
        integer :: ij

        ! solve cohort that just entered the economy
        c(1,it) = get_Psi(1, it)*get_W(1, it)
        call get_path(1,it)
        if(CES .or. CD)call shadw(1, it)

        ! derive behavior for all other cohorts in year 1 of transition
        if(it == 1)then
            do ij = 2, JJ
                c(ij,it) = get_Psi(ij,it)*get_W(ij,it)
                call get_path(ij,it)
                if(CES .or. CD)call shadw(ij,it)
            enddo

        endif

    end subroutine


    ! present value of resources for household aged ij in year it
    function get_W(ij, it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: ik, ik1, ijp
        real*8 :: get_W, PRn

        ! Compute index for year corresponding to age i
        get_W = 0.0d0
        PRn = 1.0d0

        if(it == 1) get_W=Rn(it)*a(ij, it)+v(-ij+2)
        if(ij == 1) get_W=v(it)

        do ijp = ij, JJ
            ik = year(it, ij, ijp)
            ik1 = year(it,ij, ijp+1)
            get_W = get_W +  (wn(ijp,ik) + pen(ijp,ik))/PRn
            PRn = PRn*Rn(ik1)
        enddo

    end function


    ! function for calculating the marginal consumption in the initial age
    function get_Psi(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: ik, ik1, ijp
        real*8 :: get_Psi, PRn, vv(JJ)

        get_Psi = 0d0
        PRn = 1d0

        do ijp = ij, JJ
            ik = year(it, ij, ijp)
            ik1 = year(it, ij, ijp+1)
            vv(ijp) = (1d0+nu**rho*(wn(ijp,ik)/p(ik))**(1d0-rho))**((rho-gamma)/(1d0-rho))
            if(CES)then
                get_Psi = get_Psi+beta**((ijp-ij)*gamma)*(p(ik)/PRn/p(it))**(1d0-gamma)*&
                          vv(ijp)**((1d0-gamma)/(rho-gamma))
            elseif(CD)then
                get_Psi = get_Psi+beta**((ijp-ij)*gamma)*((p(ik)/p(it))**nu/PRn)**(1d0-gamma)* &
                          (wn(ijp,ik)/wn(ij,it))**((1d0-nu)*(1d0-gamma))
            else
                get_Psi = get_Psi+beta**((ijp-ij)*gamma)*(p(ik)/PRn/p(it))**(1d0-gamma)
            endif
            PRn = PRn*Rn(ik1)
        enddo
        if(CES)then
            get_Psi = vv(ij)/p(it)/get_Psi
        elseif(CD)then
            get_Psi = nu/p(it)/get_Psi
        else
            get_Psi = 1d0/p(it)/get_Psi
        endif

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: ijp, ik, ik1, itm
        real*8 :: PRn, vv(JJ)

        if(CES)then
            l(ij, it) = (wn(ij,it)/nu/p(it))**(-rho)*c(ij,it)
        elseif(CD)then
            l(ij, it) = p(it)/wn(ij,it)*(1-nu)/nu*c(ij,it)
        endif
        vv(ij) = (1d0+nu**rho*(wn(ij,it)/p(it))**(1d0-rho))**((rho-gamma)/(1d0-rho))
        ik1 = year(it, 1, 2)
        PRn = Rn(ik1)
        do ijp = ij+1, JJ
            ik = year(it, ij, ijp)
            ik1 = year(it, ij, ijp+1)
            itm = year(it, ij, ijp-1)
            if(CES)then
                vv(ijp) = (1d0+nu**rho*(wn(ijp,ik)/p(ik))**(1d0-rho))**((rho-gamma)/(1d0-rho))
                c(ijp, ik) = vv(ijp)/vv(ij)*(beta**(ijp-ij)*PRn*p(it)/p(ik))**gamma*c(ij,it)
                l(ijp, ik) = (wn(ijp,ik)/nu/p(ik))**(-rho)*c(ijp,ik)
            elseif(CD)then
                c(ijp, ik) = (wn(ij,it)/wn(ijp,ik))**((1d0-nu)*(gamma-1d0))*(p(ik)/p(it))**(nu*(1-gamma)-1d0) &
                            *(beta**(ijp-ij)*PRn)**gamma*c(ij,it)
                l(ijp, ik) = p(ik)/wn(ijp,ik)*(1-nu)/nu*c(ijp,ik)
            else
                c(ijp, ik) = (beta**(ijp-ij)*PRn*p(it)/p(ik))**gamma*c(ij,it)
            endif
            a(ijp, ik) = wn(ijp-1,itm)*(1d0-l(ijp-1,itm)) + pen(ijp-1,itm)+ Rn(itm)*a(ijp-1,itm) - &
                        p(itm)*c(ijp-1,itm)
            if(ik == 2)a(ijp,ik) = a(ijp,ik)+v(-ijp+3)
            if(ik > 2 .and. ijp == 2) a(ijp,ik) = a(ijp,ik)+v(itm)
            PRn = PRn*Rn(ik1)
        enddo

    end subroutine


    ! compute shadow wages
    subroutine shadw(ij, it)

        implicit none

        integer, intent(in) :: ij, it
        integer :: ijp, ijp2, ik, ik2

        do ijp = ij, JJ
            ik = year(it, ij, ijp)
            if(l(ijp, ik) > 1d0 - 5d-2)then
                do ijp2 = ijp, JJ
                    ik2 = year(ik, ijp, ijp2)
                    if(CES)then
                        mu(ijp2, ik2) = (1-damp)*mu(ijp2, ik2)+damp*(c(ijp2,ik2)**(1d0/rho) &
                                *nu*p(ik2)/(1d0-tauw(ik2)-taup(ik2)) - h(ijp2)*w(ik2))
                    elseif(CD)then
                        mu(ijp2, ik2) = (1-damp)*mu(ijp2, ik2)+damp*((p(ik2)*(1-nu)/nu*c(ijp2,ik2)) &
                                /(1d0-tauw(ik2)-taup(ik2)) - h(ijp2)*w(ik2))
                    endif
                    if(mu(ijp2, ik2) < 0) mu(ijp2, ik2) = 0
                enddo
                exit
            endif
            mu(ijp, ik) = (1-damp)*mu(ijp, ik)
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
        integer :: it, it1, ijp
        real*8 :: sum1, sum2, ustar, PVV

        ! calculate utility for each generation
        do it = 1, TT
            call utility(it)
        enddo

        ! transfers to old generations
        BA = 0
        do ijp = 2, JJ
            v(-ijp+2) = v(-ijp+2)+get_W(1,ijp)*((util(ijp,0)/util(ijp,1))**(1d0/egam)-1d0)
            BA(2) = BA(2)+v(-ijp+2)*m(ijp,1)
        enddo

        ! long run equilibrium
        PVV  = v(TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum1 = get_W(1,TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum2 = get_W(1,TT)*(util(1,TT)*egam)**(-1d0/egam)*(1d0+r(TT))/ &
                (r(TT)-n_p(TT))

        ! transition path
        do it = TT-1, 1, -1
            it1 = year(it, 1, 2)
            PVV = PVV*(1d0+n_p(it1))/(1d0+r(it1)) + v(it)
            sum1 = sum1*(1d0+n_p(it1))/(1d0+r(it1)) + get_W(1,it)
            sum2 = sum2*(1d0+n_p(it1))/(1d0+r(it1)) + get_W(1,it)* &
                   (util(1,it)*egam)**(-1d0/egam)
        enddo

        ! calculate ustar for future generations
        ustar = ((sum1-BA(2)-PVV)/sum2)**egam/egam

        ! calculate transfers to future generations and debt of LSRA
        do it = 1, TT
            v(it) = v(it) + get_W(1, it)*((ustar/util(1, it))**(1d0/egam)-1d0)
            if(it == 2) BA(2) = (BA(2) + v(1))/(1d0+n_p(2))
            if(it > 2) BA(it) = ((1d0+r(it-1))*BA(it-1) + v(it-1))/(1d0+n_p(it))
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, it1, ijp

        itm = year(it, 2, 1)
        it1 = year(it, 1, 2)

        if(it == 0)then
            GG(it) = gy*YY(it)
        else
            GG(it) = GG(0)
        endif

        CC(it) = 0d0
        AA(it) = 0d0
        LL(it) = 0d0
        do ijp = 1, JJ
            CC(it) = CC(it)+c(ijp,it)*m(ijp,it)
            AA(it) = AA(it)+a(ijp,it)*m(ijp,it)
            LL(it) = LL(it)+h(ijp)*(1d0 - l(ijp, it))*m(ijp,it)
        enddo
        YY(it) = KK(it)**alpha * LL(it)**(1d0-alpha)
        BB(it) = by(itm)*YY(it)

        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)*(1d0-eps(it)*tauk(it))/(1d0-tauk(it))+ &
                    delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BA(it) - BB(it)
            TB(it) = (1d0+n_p(it1))*BF(it1) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)-BA(it)) + (1d0-damp)*KK(it)
        endif
        II(it) = (1d0+n_p(it1))*KK(it1)-(1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: it1, ijp
        real*8 :: taxrev(4)

        it1 = year(it, 1, 2)

        taxrev(1) = tauc(it)*CC(it)
        taxrev(2) = tauw(it)*w(it)*LL(it)
        taxrev(3) = taur(it)*r(it)*AA(it)
        taxrev(4) = tauk(it)*(YY(it)-w(it)*LL(it)-(delta+eps(it)*r(it))*KK(it))

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(2) + taxrev(3) + &
                        taxrev(4) + (1d0+n_p(it1))*BB(it1)))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(4)+&
                        (1d0+n_p(it1))*BB(it1)))/(w(it)*LL(it)+r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(3) + &
                        taxrev(4) + (1d0+n_p(it1))*BB(it1)))/(w(it)*LL(it))
        else
            taur(it) = ((1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(2) + &
                        taxrev(4) + (1d0+n_p(it1))*BB(it1)))/(r(it)*AA(it))
        endif
        TXR(it) = sum(taxrev)
        pen(JR:JJ,it) = kappa(it)*w(it)

        ! get budget balancing social security contribution
        Tpen(it) = 0d0
        do ijp = JR, JJ
            Tpen(it) = Tpen(it) + pen(ijp,it)*M(ijp,it)
        enddo
        taup(it) = Tpen(it)/w(it)/LL(it)

    end subroutine


    ! subroutine for writing output to a file or the console (fname = 0)
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff
        integer :: it1, ijp

        ! household utility
        call utility(it)

        it1 = year(it, 1, 2)

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
        write(fname,'(3f7.2/)')LL(it),w(it),util(1,it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk    TXR     DD     rB  '
        write(fname,'(8f7.2/)')tauc(it),tauw(it),taur(it),taup(it),tauk(it),TXR(it)/YY(it),&
                            ((1+n_p(it1))*BB(it1)-BB(it))/YY(it),r(it)*BB(it)/YY(it)


        if(fname > 1)then
            write(fname,'(a)')' Age    cons     leis     wn       mu       pen      asset    Diff      '
            do ijp = 1, JJ
                it1 = year(it, 1, 2)
                if(ijp < JJ)then
                    diff = Rn(it)*a(ijp,it)+wn(ijp,it)*(1d0-l(ijp,it))+pen(ijp,it)-a(ijp+1,it1)-p(it)*c(ijp,it)
                else
                    diff = Rn(it)*a(ijp,it)+wn(ijp,it)*(1d0-l(ijp,it))+pen(ijp,it)-p(it)*c(ijp,it)
                endif
                write(fname,'(i3,7f9.2)')ijp,c(ijp,it),l(ijp,it),wn(ijp,it),mu(ijp,it),pen(ijp,it),a(ijp,it),diff
            enddo
            write(fname,'(a/)')'                                           '
        endif

    end subroutine


    ! subroutine for writing output summary to a file or the console (fname = 0)
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it, ijp
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'         C      A      K      L      r', &
            '      w      Y      B     BA   tauc   tauw   taur', &
            '   taup   tauk      v    HEV   y-d'
        diff = YY(0)-CC(0)-GG(0)-II(0)
        ! current generations
        do ijp = JJ, 3, -1
            write(fname,'(i3,98x,2f7.2)')-ijp+2, v(-ijp+2), &
                ((util(ijp,1)/util(ijp,0))**(1d0/egam)-1d0)*100d0
        enddo
        write(fname,'(i3,17f7.2)')0,CC(0),AA(0),KK(0),LL(0),r(0),w(0), &
            YY(0),BB(0),BA(0),tauc(0),tauw(0),taur(0),taup(0),tauk(0),v(0), &
            ((util(2,1)/util(2,0))**(1d0/egam)-1d0)*100d0,diff

        ! future generations
        do it = 1, TT
            diff = YY(it)-CC(it)-GG(it)-II(it)
            write(fname,'(i3,17f7.2)')it,CC(it),AA(it),KK(it),LL(it),r(it),w(it), &
                YY(it),BB(it),BA(it),tauc(it),tauw(it),taur(it),taup(it),tauk(it),v(it), &
                ((util(1,it)/util(1,0))**(1d0/egam)-1d0)*100d0,diff
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  ik, itp, ijp

        util(1,it) = 0d0
        do ijp = 1, JJ
            ik = year(it, 1, ijp)
            if(CES)then
                util(1,it) = util(1,it)+beta**(ijp-1)*(c(ijp,ik)**erho+nu*l(ijp,ik)**erho)**(egam/erho)/egam
            elseif(CD)then
                util(1,it) = util(1,it)+beta**(ijp-1)*(c(ijp,ik)**nu*l(ijp,ik)**(1d0-nu))**egam/egam
            else
                util(1,it) = util(1,it)+beta**(ijp-1)*c(ijp,ik)**egam/egam
            endif
        enddo

        if(it < 2)then
            do itp = 2, JJ
                util(itp,it) = 0d0
                do ijp = itp, JJ
                    ik = year(it, itp, ijp)
                    if(CES)then
                        util(itp,it) = util(itp,it) + beta**(ijp-itp)*(c(ijp,ik)**erho + &
                                               nu*l(ijp,ik)**erho)**(egam/erho)/egam
                    elseif(CD)then
                        util(itp,it) = util(itp,it) + beta**(ijp-itp)*(c(ijp,ik)**nu*l(ijp,ik)**(1d0-nu)) &
                                                                        **egam/egam
                    else
                        util(itp,it) = util(itp,it) + beta**(ijp-itp)*c(ijp,ik)**egam/egam
                    endif
                enddo
            enddo
        endif

    end subroutine

end program
