!##############################################################################
! PROGRAM ENGROW_OLG
!
! ##  Endogenous growth economy
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program ENGROW_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 36    ! 160
    integer, parameter :: JJ      = 6       !  55
    integer, parameter :: JR      = 5       ! 41
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  delh    = 0.5d0
    real*8, parameter ::  xi      = 2.0d0
    real*8, parameter ::  upsi    = 0.5d0
    real*8, parameter ::  epsi    = 0.0d0
    real*8, parameter ::  mu      = 1.0d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.4d0
    integer, parameter :: itermax = 100

    ! model variables
    real*8 :: w(0:TT), r(0:TT), wn(0:TT), winc(JJ,0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT), pen(JJ,0:TT)
    real*8 :: by(0:TT), kappa(0:TT), n_p(0:TT), gy, tauk(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BF(0:TT), TB(0:TT)
    real*8 :: Tpen(0:TT), TXR(0:TT), m(JJ,0:TT), ma(JJ,0:TT)
    real*8 :: a(JJ,0:TT), c(JJ,0:TT), util(JJ,0:TT)
    real*8 :: taus(0:TT), e(0:TT), HH(0:TT), ne(0:TT)
    integer :: tax(0:TT), eps(0:TT)
    logical :: smopec
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


contains

    ! initializes variables and government parameters
    subroutine initialize

        implicit none
        integer :: ij, it, itm

        ! set model parameters
        n_p(0) = 0d0
        n_p(1:TT) = 0d0
        gy = 0.0d0
        smopec = .false.

        ! set reform values
        by(0:TT) = 0d0
        kappa(0) = 0.0d0
        kappa(1:TT) = 0.0d0
        eps(0) = 0
        eps(1:TT) = 0
        tauk(0) = 0d0
        tauk(1:TT) = 0d0
        taus(0) = 0d0
        taus(1:TT) = 0d0
        tax(0) = 1
        tax(1:TT) = 1

        ! initialize tax rates and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0

        ! initialize assets and debt holdings
        a = 0d0
        e = 0d0
        wn = 0d0
        winc = 0d0
        ne = 0d0
        YY = 0d0
        HH = 1d0
        BF = 0d0
        TB = 0d0
        TXR = 0d0
        ma = 0d0


        ! size of cohorts in specific year
        do it = 0, TT
            m(1,it) = 1d0
            itm = year(it, 2, 1)
            do ij = 2, JJ
                m(ij,it) = m(ij-1,itm)/(1+n_p(it))
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
        a(:,1) = a(:,0)
        ne(1) = ne(0)

        KK(:) = KK(0)
        LL(:) = LL(0)
        HH(:) = HH(0)


        do iter = 1, itermax

            ! get prices, decisions and quantities
            do it = 1, TT
                call factor_prices(it)
            enddo

            do it = 1, TT
                call decisions(it)
            enddo

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
        integer :: itp, ij

        ! factor prices and pension payments in year t
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = (1d0-tauk(it))/(1d0-eps(it)*tauk(it))*     &
                    (alpha*(KK(it)/LL(it))**(alpha-1d0)*HH(it)**epsi-delta)
        endif
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha*HH(it)**epsi
        wn(it) = w(it)*(1d0-tauw(it)-taup(it))
        winc(1,it) = (1d0-(1d0-taus(it))*e(it))*wn(it)
        do ij = 2, JR-1
            itp = year(it, ij, 2)
            winc(ij,it) = wn(it)*(1d0-delh)**(ij-2)*(1+ne(itp))/mu
        enddo
        Rn(it) = 1d0+r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it
        integer :: ij


        call getedu(it)

        c(1,it) = get_Psi(1,it)*get_W(1,it)
        call get_path(1,it)

        if(it == 1)then
            do ij = 2, JJ
                c(ij,it) = get_Psi(ij,it)*get_W(ij,it)
                call get_path(ij,it)
            enddo
        endif

    end subroutine


    ! subroutine for calculation education investment and human capital
    subroutine getedu(it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it
        integer :: ij, itp
        real*8 :: pvw, PRn


        ! Compute pv of net wage
        pvw = 0d0
        PRn = 1d0
        do ij = 2, JR-1
            itp = year(it, 1, ij)
            PRn = PRn*Rn(itp)
            pvw = pvw + wn(itp)*(1d0-delh)**(ij-2)/PRn
        enddo

        e(it) = (xi*upsi*pvw/(wn(it)*(1d0-taus(it))))**(1d0/(1d0-upsi))
        if(e(it) > 0.99d0)e(it) = 0.99d0

        itp = year(it, 1, 2)
        ne(itp) = mu*(1d0+xi*e(it)**upsi)-1d0

    end subroutine


    ! function for calculation present value of resources for household age j in year it
    function get_W(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: itp, ijp
        real*8 :: get_W, PRn

        ! Compute index for year corresponding to age i
        get_W = winc(ij,it)+pen(ij,it)
        if(it == 1 .and. ij > 1) get_W = get_W + Rn(it)*a(ij,it)

        PRn = 1d0
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_W = get_W +  (winc(ijp,itp)+pen(ijp,itp))/PRn
        enddo

    end function


    ! function for calculating the marginal consumption in the initial age
    function get_Psi(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: itp, ijp
        real*8 :: get_Psi, PRn

        get_Psi = 1d0
        PRn = 1d0
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_Psi = get_Psi+beta**((ijp-ij)*gamma)*(p(itp)/PRn/p(it))**(1d0-gamma)
        enddo
        get_Psi = 1d0/get_Psi/p(it)

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: ijp, itp, itm
        real*8 :: PRn

        PRn = 1d0
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            itm = year(it, ij, ijp-1)
            c(ijp, itp) = (beta**(ijp-ij)*PRn*p(it)/p(itp))**gamma*c(ij,it)
            a(ijp, itp) = winc(ijp-1,itm) + pen(ijp-1,itm)+ Rn(itm)*a(ijp-1,itm) - p(itm)*c(ijp-1,itm)
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


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp, ij

        itm = year(it, 2, 1)

        if(it == 0)then
            GG(it) = gy*YY(it)
        else
            GG(it) = GG(0)
        endif

        ma(1,it) = 1d0
        do ij = 2, JJ
            ma(ij,it) =  ma(ij-1,itm)/(1d0+ne(it))/(1d0+n_p(it))
        enddo


        CC(it) = 0d0
        AA(it) = 0d0
        LL(it) = 1d0-e(it)
        HH(it) = 0d0
        do ij = 1, JJ
            itp = year(it, 2, ij)
            CC(it) = CC(it)+c(ij,it)*ma(ij,it)
            AA(it) = AA(it)+a(ij,it)*ma(ij,it)
            if(ij > 1 .and. ij < JR) LL(it) = LL(it)+(1-delh)**(ij-2)*(1d0+ne(itp))/mu*ma(ij,it)
            if(ij < JR) HH(it) = HH(it)+m(ij,it)
        enddo
        HH(it) = LL(it)/(HH(it) - e(it))

        YY(it) = KK(it)**alpha*LL(it)**(1d0-alpha)*HH(it)**epsi
        BB(it) = by(itm)*YY(it)


        itp = year(it, 1, 2)
        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)*(1d0-eps(it)*tauk(it))/(1d0-tauk(it))+ &
                    delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BB(it)
            TB(it) = (1d0+n_p(itp))*(1d0+ne(itp))*BF(itp) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)) + (1d0-damp)*KK(it)
        endif
        II(it) = (1d0+n_p(itp))*(1d0+ne(itp))*KK(itp)-(1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: it1, ij
        real*8 :: taxrev(5)

        it1 = year(it, 1, 2)

        taxrev(1) = tauc(it)*CC(it)
        taxrev(2) = tauw(it)*w(it)*LL(it)
        taxrev(3) = taur(it)*r(it)*AA(it)
        taxrev(4) = tauk(it)*(YY(it)-w(it)*LL(it)-(delta+eps(it)*r(it))*KK(it))
        taxrev(5) = taus(it)*e(it)*wn(it)

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = (taxrev(5)+(1d0+r(it))*BB(it)+GG(it)-(taxrev(2) + taxrev(3)+&
                        taxrev(4)+(1d0+n_p(it1))*(1d0+ne(it1))*BB(it1)))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = (taxrev(5)+(1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(4)+&
                        (1d0+n_p(it1))*(1d0+ne(it1))*BB(it1)))/(w(it)*LL(it)+r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = (taxrev(5)+(1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(3) + &
                        taxrev(4) + (1d0+n_p(it1))*(1d0+ne(it1))*BB(it1)))/(w(it)*LL(it))
        else
            taur(it) = (taxrev(5)+(1d0+r(it))*BB(it)+GG(it)-(taxrev(1) + taxrev(2) + &
                        taxrev(4) + (1d0+n_p(it1))*(1d0+ne(it1))*BB(it1)))/(r(it)*AA(it))
        endif
        TXR(it) = sum(taxrev(1:4))-taxrev(5)


        pen(JR:JJ,it) = kappa(it)*winc(JR-1,it)
        ! get budget balancing social security contribution
        Tpen(it) = 0d0
        do ij = JR , JJ
            Tpen(it) = Tpen(it) + pen(ij,it)*ma(ij,it)
        enddo
        taup(it) = Tpen(it)/w(it)/LL(it)

    end subroutine


    ! subroutine for writing output to a file or the console (fname = 0)
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff
        integer :: it1, ij

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
        write(fname,'(a)')'      A      K     BB     BF       r'
        write(fname,'(4f7.2,f8.2)')AA(it), KK(it), BB(it), BF(it), r(it)
        write(fname,'(4f7.2/)')AA(it)/YY(it),KK(it)/YY(it),BB(it)/YY(it),BF(it)/YY(it)
        write(fname,'(a)')' Labor Market  '
        write(fname,'(a)')'     LL      e     ne      w   util'
        write(fname,'(5f7.2/)')LL(it),e(it),ne(it),w(it),util(1,it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk   taus    TXR     DD     rB  '
        write(fname,'(9f7.2/)')tauc(it),tauw(it),taur(it),taup(it),tauk(it),taus(it),TXR(it)/YY(it),&
                           ((1+n_p(it1))*BB(it1)-BB(it))/YY(it),r(it)*BB(it)/YY(it)


        if(fname > 1)then
            write(fname,'(a)')' Age    cons       wn      pen    asset    Diff      '
            do ij = 1, JJ
                it1 = year(it, 1, 2)
                if(ij < JJ)then
                    diff = Rn(it)*a(ij,it)+winc(ij,it)+pen(ij,it)-a(ij+1,it1)-p(it)*c(ij,it)
!                   if(ij == 1) diff = diff - (1d0-taus(it))*e(it)*wn(it)
                    write(fname,'(i3,5f9.2)')ij,c(ij,it),winc(ij,it),pen(ij,it),a(ij,it),diff
                else
                    diff = Rn(it)*a(ij,it)+pen(ij,it)-p(it)*c(ij,it)
                    write(fname,'(i3,5f9.2)')ij,c(ij,it),0.00,pen(ij,it),a(ij,it),diff
                endif
            enddo
            write(fname,'(a/)')'                                           '
        endif

    end subroutine


    ! subroutine for writing output summary to a file or the console (fname = 0)
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it, ij
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'         C      e     ne      L      K      r', &
            '      w      Y      B   tauc   tauw   taur', &
            '   taup   tauk   taus    HEV   y-d'
        diff = YY(0)-CC(0)-GG(0)-II(0)
        ! current generations
        do ij = JJ, 3, -1
            write(fname,'(i3,112x,f7.2)')-ij+2, &
                    ((util(ij,1)/util(ij,0))**(1d0/egam)-1d0)*100d0
        enddo
        write(fname,'(i3,17f7.2)')0,CC(0),e(0),ne(0),LL(0),KK(0),r(0),w(0), &
            YY(0),BB(0),tauc(0),tauw(0),taur(0),taup(0),tauk(0),taus(0),&
            ((util(2,1)/util(2,0))**(1d0/egam)-1d0)*100d0,diff

        ! future generations
        do it = 1, TT
            diff = YY(it)-CC(it)-GG(it)-II(it)
            write(fname,'(i3,17f7.2)')it,CC(it),e(it),ne(it),LL(it),KK(it),r(it),w(it), &
                    YY(it),BB(it),tauc(it),tauw(it),taur(it),taup(it),tauk(it),taus(it),&
                    ((util(1,it)/util(1,0))**(1d0/egam)-1d0)*100d0,diff
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  ik, iu, ij

        util(1,it) = 0d0
        do ij = 1, JJ
            ik = year(it, 1, ij)
            util(1,it) = util(1,it)+beta**(ij-1)*c(ij,ik)**egam/egam
        enddo

        if(it < 2)then
            do iu = 2, JJ
                util(iu,it) = 0d0
                do ij = iu, JJ
                    ik = year(it, iu, ij)
                    util(iu,it) = util(iu,it) + beta**(ij-iu)*c(ij,ik)**egam/egam
                enddo
            enddo
        endif

    end subroutine

end program
