!##############################################################################
! PROGRAM TRHCEG_OLG
!
! ## The OLG model with endogenous growth
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
program TRHCEG_OLG

    implicit none

    ! model parameters
    integer, parameter :: TT      = 24
    integer, parameter :: JJ      = 3
    integer, parameter :: JR      = 3
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  xi      = 2.0d0
    real*8, parameter ::  upsi    = 0.5d0
    real*8, parameter ::  epsi    = 0.0d0
    real*8, parameter ::  mu      = 0.8d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.3d0
    integer, parameter :: itermax = 200

    ! model variables
    real*8 :: w(0:TT), r(0:TT), wn(0:TT), winc(JJ,0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT), pen(JJ,0:TT)
    real*8 :: by(0:TT), kappa(0:TT), n_p(0:TT), gy, tauk(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BF(0:TT), TB(0:TT)
    real*8 :: Tpen(0:TT), TXR(0:TT), m(JJ,0:TT), ma(JJ,0:TT)
    real*8 :: a(JJ,0:TT), c(JJ,0:TT)
    real*8 :: taus(0:TT), e(0:TT), HH(0:TT), n_e(0:TT)
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
        n_p = 0d0
        by = 0d0
        kappa = 0d0
        eps = 0
        tauk = 0d0
        taus = 0d0
        tax = 1
        smopec = .false.


        ! impact of mu Table 7.8
        ! change mu in the parameter section at beginning of program
        gy = 0.00d0


        ! setup for Table 7.9
        ! set mu = 0.8d0 at beginnig of program section
        !gy = 0.22d0
        !tax(1:TT) = 3


        ! setup for Table 7.9
        ! set mu = 0.8d0 at beginnig of program section
        !gy = 0.22d0
        !taus(1:TT) = 0.25d0


        ! initialize tax rates and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0

        ! initialize assets, LSRA payments and debt holdings
        a = 0d0
        e = 0d0
        winc = 0d0
        n_e = 0d0
        YY = 0d0
        HH = 1d0
        BF = 0d0
        TB = 0d0
        TXR = 0d0
        ma = 0d0

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
        n_e(1) = n_e(0)
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
        integer :: ij, itp

        ! factor prices and pension payments in year t
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = (1d0-tauk(it))/(1d0-eps(it)*tauk(it))*     &
                    (alpha*(KK(it)/LL(it))**(alpha-1d0)*HH(it)**epsi-delta)
        endif
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha*HH(it)**epsi

        ! derive after tax prices
        wn(it) = w(it)*(1d0-tauw(it)-taup(it))
        winc(1, it) = (1d0 - (1d0-taus(it))*e(it))*wn(it)
        do ij = 2, JR-1
            itp = year(it, ij, 2)
            winc(ij, it) = wn(it)*(1d0+n_e(itp))/mu
        enddo
        Rn(it) = 1d0 + r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none
        integer, intent(in) :: it
        integer :: ij, itp

        ! derive education decision and corresponding growth rate
        itp = year(it, 1, 2)
        e(it) = (xi*upsi*wn(itp)/(wn(it)*(1d0-taus(it))*Rn(itp)))**(1d0/(1d0-upsi))
        n_e(itp) = mu*(1d0+xi*e(it)**upsi) - 1d0

        ! consumption path for cohort that just entered the economy
        c(1, it) = get_Psi(1, it)*get_W(1, it)
        call get_path(1, it)

        ! consumption for all other cohorts in year 1 of transition
        if(it == 1)then
            do ij = 2, JJ
                c(ij, it) = get_Psi(ij, it)*get_W(ij, it)
                call get_path(ij, it)
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
        get_W = winc(ij, it) + pen(ij, it)
        if(it == 1 .and. ij > 1)get_W = get_W + Rn(it)*a(ij, it)

        ! iterate over remainder of life span
        PRn = 1d0
        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_W = get_W + (winc(ijp, itp) + pen(ijp, itp))/PRn
        enddo

    end function


    ! marginal consumption for household aged ij in year it
    function get_Psi(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: itp, ijp
        real*8 :: get_Psi, PRn

        get_Psi = 1d0
        PRn = 1d0

        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_Psi = get_Psi + beta**((ijp-ij)*gamma)*(p(itp)/PRn/p(it))**(1d0-gamma)
        enddo
        get_Psi = 1d0/p(it)/get_Psi

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ij, it)

        implicit none
        integer, intent(in) :: ij, it
        integer :: ijp, itp, itm
        real*8 :: PRn

        PRn = 1d0

        ! determine consumption and leisure path for remainder of lifetime
        do ijp = ij+1, JJ

            ! get future and previous year as well as interest factor
            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp-1)
            PRn = PRn*Rn(itp)

            ! get consumption and assets
            c(ijp, itp) = (beta**(ijp-ij)*PRn*p(it)/p(itp))**gamma*c(ij, it)
            a(ijp, itp) = winc(ijp-1, itm) + pen(ijp-1, itm) + &
                          Rn(itm)*a(ijp-1, itm) - p(itm)*c(ijp-1, itm)
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp, ij

        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        ! government consumption
        if(it == 0)then
            GG(it) = gy*YY(it)
        else
            GG(it) = GG(0)
        endif

        ! calculate adjusted aggregator
        ma(1, it) = 1d0
        do ij = 2, JJ
            ma(ij, it) = ma(ij-1, itm)/((1d0+n_e(it))*(1d0+n_p(it)))
        enddo

        ! calculate economy wide aggregates
        CC(it) = 0d0
        AA(it) = 0d0
        LL(it) = 0d0
        HH(it) = 0d0
        do ij = 1, JJ
            CC(it) = CC(it) + c(ij, it)*ma(ij, it)
            AA(it) = AA(it) + a(ij, it)*ma(ij, it)
            if(ij < JR)LL(it) = LL(it) + mu**dble(1-ij)*m(ij, it)
            if(ij < JR)HH(it) = HH(it) + m(ij, it)
        enddo
        LL(it) = LL(it) - e(it)
        HH(it) = LL(it)/(HH(it) - e(it))

        ! derive output and government debt
        YY(it) = KK(it)**alpha*LL(it)**(1d0-alpha)*HH(it)**epsi
        BB(it) = by(itm)*YY(it)

        ! capital stock (in closed or open economy)
        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)*(1d0-eps(it)*tauk(it))/(1d0-tauk(it))+ &
                         delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BB(it)
            TB(it) = (1d0+n_p(itp))*(1d0+n_e(itp))*BF(itp) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)) + (1d0-damp)*KK(it)
        endif
        II(it) = (1d0+n_p(itp))*(1d0+n_e(itp))*KK(itp)-(1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: itp, ij
        real*8 :: taxrev(5), n_g

        itp = year(it, 1, 2)

        taxrev(1) = tauc(it)*CC(it)
        taxrev(2) = tauw(it)*w(it)*LL(it)
        taxrev(3) = taur(it)*r(it)*AA(it)
        taxrev(4) = tauk(it)*(YY(it) - w(it)*LL(it) - (delta+eps(it)*r(it))*KK(it))
        taxrev(5) = taus(it)*e(it)*wn(it)

        ! total growth rate
        n_g = (1d0+n_p(itp))*(1d0+n_e(itp)) - 1d0

        ! get budget balancing tax rate
                if(tax(it) == 1)then
            tauc(it) = (taxrev(5) + (1d0+r(it))*BB(it) + GG(it) - (taxrev(2) + &
                        taxrev(3) + taxrev(4) + (1d0+n_g)*BB(itp)))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = (taxrev(5) + (1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + &
                        taxrev(4) + (1d0+n_g)*BB(itp)))/(w(it)*LL(it) + r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = (taxrev(5) + (1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + &
                        taxrev(3) + taxrev(4) + (1d0+n_g)*BB(itp)))/(w(it)*LL(it))
        else
            taur(it) = (taxrev(5) + (1d0+r(it))*BB(it) + GG(it) - (taxrev(1) + &
                        taxrev(2) + taxrev(4) + (1d0+n_g)*BB(itp)))/(r(it)*AA(it))
        endif
        TXR(it) = sum(taxrev(1:4)) - taxrev(5)

        ! get budget balancing social security contribution
        pen(JR:JJ, it) = kappa(it)*w(it)
        Tpen(it) = 0d0
        do ij = JR, JJ
            Tpen(it) = Tpen(it) + pen(ij, it)*ma(ij, it)
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
        real*8 :: diff, n_g
        integer :: itp, ij

        itp = year(it, 1, 2)

        ! total growth rate
        n_g = (1d0+n_p(itp))*(1d0+n_e(itp)) - 1d0

        ! difference on goods market
        diff = YY(it) - CC(it) - GG(it) - II(it) - TB(it)

        write(fname,'(a,i3/)')'Equilibrium: Year', it
        write(fname,'(a)')' Goods Market  '
        write(fname,'(a)')'      Y      C      G      I     TB      DIFF'
        write(fname,'(5f7.2,f10.4)')YY(it), CC(it), GG(it), II(it), TB(it), diff
        write(fname,'(5f7.2/)')(/YY(it), CC(it), GG(it), II(it), TB(it)/)/YY(it)*100d0
        write(fname,'(a)')' Capital Market  '
        write(fname,'(a)')'      A      K     BB     BF       r'
        write(fname,'(4f7.2,f8.2)')AA(it), KK(it), BB(it), BF(it), r(it)
        write(fname,'(4f7.2/)')(/AA(it), KK(it), BB(it), BF(it)/)/YY(it)*100d0
        write(fname,'(a)')' Labor Market  '
        write(fname,'(a)')'     LL      e    n_e      w'
        write(fname,'(4f7.2/)')LL(it), e(it), n_e(it), w(it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk   taus    TXR     DD     rB  '
        write(fname,'(9f7.2/)')tauc(it), tauw(it), taur(it), taup(it), tauk(it), taus(it), &
                        (/TXR(it), (1d0+n_g)*BB(itp)-BB(it), r(it)*BB(it)/)/YY(it)*100d0

        ! write individual output
        write(fname,'(a)')' Age    cons       wn      pen    asset    Diff      '
        do ij = 1, JJ
            itp = year(it, 1, 2)
            if(ij < JJ)then
                diff = Rn(it)*a(ij, it) + wn(it)*((1d0+n_e(it))/mu)**(ij-1) &
                       + pen(ij, it) - a(ij+1, itp) - p(it)*c(ij, it)
                if(ij == 1)diff = diff - (1d0-taus(it))*e(it)*wn(it)
                write(fname,'(i3,5f9.2)')ij, c(ij, it), wn(it)*((1d0+n_e(it))/mu)**(ij-1),&
                                         pen(ij, it), a(ij, it), diff
            else
                diff = Rn(it)*a(ij, it) + pen(ij, it) - p(it)*c(ij, it)
                write(fname,'(i3,5f9.2)')ij, c(ij, it), 0d0, pen(ij, it), a(ij, it), diff
            endif
        enddo
        write(fname,'(/)')

    end subroutine


    ! subroutine for writing output summary to a file
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'         C      e    n_e      L      K      r', &
            '      w      Y      B   tauc   tauw   taur', &
            '   taup   tauk   taus   y-d'

        diff = YY(0) - CC(0) - GG(0) - II(0)

        ! current generations
        write(fname,'(i3,16f7.2)')0, CC(0), e(0), n_e(0), LL(0), KK(0), r(0), w(0), &
            YY(0), BB(0), tauc(0), tauw(0), taur(0), taup(0), tauk(0), taus(0), diff

        ! future generations
        do it = 1, TT
            diff = YY(it) - CC(it) - GG(it) - II(it)
            write(fname,'(i3,16f7.2)')it, CC(it), e(it), n_e(it), LL(it), KK(it), r(it), w(it), &
                YY(it), BB(it), tauc(it), tauw(it), taur(it), taup(it), tauk(it), taus(it), diff
        enddo

    end subroutine

end program
