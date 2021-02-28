!##############################################################################
! PROGRAM ULTSIXBEQ_OLG
!
! ##  uncertain lifetime, six periods and bequest
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
include "sol_prog07_11m.f90"

program ULTSIXBEQ_OLG

    use globals
    use toolbox

    implicit none

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
        integer :: ij, it, itm

        ! set model parameters
        n_p(0) = 0.2d0
        n_p(1:TT) = 0.2d0
        gy = 0.20d0
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

        ! initialize tax rates and pensions
        tauc = 0d0
        tauw = 0d0
        taur = 0d0
        taup = 0d0
        pen = 0d0

        ! initialize assets, LSRA payments and debt holdings
        a = 0d0
        v = 0d0
        beq = 0d0
        BQ = 0d0
        GAM = 0d0
        YY = 0d0
        BA = 0d0
        BF = 0d0
        TB = 0d0
        TXR = 0d0
        util = 0d0

        ! human capital profile
        do ij = 1, JR-1
!           h(ij) = 1.0d0 + 0.05*ij-0.001*ij**2
            h(ij) = 1d0
        enddo
        h(JR:JJ) = 0d0

        ! set survival probs
        psi(1,:) = 1.0d0
        psi(2,0:TT) = 0.95d0
        psi(3,0:TT) = 0.85d0
        psi(4,0:TT) = 0.75d0
        psi(5,0:TT) = 0.7d0
        psi(6,0:TT) = 0.65d0
        psi(7,0:TT) = 0.0d0
!       psi(1,:) = 1.0d0
!       psi(2,0:TT) = 0.95d0
!       psi(3,0:TT) = 0.95d0
!       psi(4,0:TT) = 0.90d0
!       psi(5,0:TT) = 0.85d0
!       psi(6,0:TT) = 0.80d0
!       psi(7,0:TT) = 0.0d0

        ! set bequest distribution
        omebeq(1) = 0.5d0
        omebeq(2) = 0.5d0
        omebeq(3) = 0.0d0
        omebeq(4) = 0.0d0
        omebeq(5) = 0d0
        omebeq(6) = 0d0


        ! size of cohorts in specific year
        do it = 0, TT
            m(1,it) = 1d0
            GAM(1,it) = omebeq(1)
            itm = year(it, 2, 1)

            do ij = 2, JJ
                m(ij,it) = m(ij-1,itm)*psi(ij,it)/(1+n_p(it))
                GAM(1,it) = GAM(1,it) + omebeq(ij)*m(ij,it)
            enddo
            do ij = JJ, 1, -1
                GAM(ij,it) = omebeq(ij)/GAM(1,it)
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

        write(*,'(a,i4,a,i4,a,f16.10)')'Iteration: ',iter,' Markets: ', nmarket,' Diff: ', &
                                   maxval(abs(YY(:) - CC(:) - II(:) - GG(:) - TB(:))/YY(:))
        write(*,*)
        if(iter > itermax)then
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        endif

    end subroutine


    ! subroutine for calculating factor prices in a certain year
    subroutine factor_prices(it)

        implicit none
        integer, intent(in) :: it

        ! factor prices and pension payments in year t
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = (1d0-tauk(it))/(1d0-eps(it)*tauk(it))*     &
                        (alpha*(KK(it)/LL(it))**(alpha-1d0)-delta)
        endif
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha
        wn(it) = w(it)*(1d0-tauw(it)-taup(it))
        Rn(it) = 1d0+r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it
        integer :: ij

        call get_path(1,it)

        if(it == 1)then
            do ij = 2, JJ
                call get_path(ij,it)
            enddo
        endif

    end subroutine


    ! function for calculation present value of resources for household age j in year it
    function get_W(ij,it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: itp, ijp
        real*8 :: get_W, PRn

        ! Compute index for year corresponding to age i
        get_W = wn(it)*h(ij) + beq(ij,it) + pen(ij,it)
        PRn = 1.0d0

        if(it == 1 .and. ij > 1) get_W = get_W + Rn(it)*a(ij,it)+v(2-ij)
        if(it >= 1 .and. ij == 1)get_W = get_W + v(it)

        do ijp = ij+1, JJ
            itp = year(it, ij, ijp)
            PRn = PRn*Rn(itp)
            get_W = get_W +  (wn(itp)*h(ijp) + beq(ijp,itp) + pen(ijp,itp))/PRn
        enddo

    end function


    ! subroutine for calculating the optimal consumption path
    subroutine get_path(ij, it)

        implicit none

        ! variable declaration
        integer, intent(in) :: it, ij
        integer :: ijp, itp, itm
        real*8 :: x(JJ-ij+1), vh
        logical :: check

        it_com = it
        ij_com = ij

        x = 0.2d0
        call fzero(x, foc, check)

        if(check) stop 'fzero did not converge'

        do ijp = ij, JJ

            vh = 0d0
            if(it == 1 .and. ij > 1)vh = v(2-ij)
            if(it >= 1 .and. ij == 1)vh = v(it)

            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp+1)
            c(ijp, itp) = x(ijp-ij+1)
            a(ijp+1, itm) = Rn(itp)*a(ijp,itp) + wn(itp)*h(ijp) + &
                beq(ijp,itp) + pen(ijp,itp) + vh - p(itp)*c(ijp, itp)
        enddo

    end subroutine


    ! subroutine for calculating lsra transfers
    subroutine lsra()

        implicit none
        integer :: it, it1, ij
        real*8 :: sum1, sum2, ustar, PVV

        ! calculate utility for each generation
        do it = 1, TT
            call utility(it)
        enddo

        ! transfers to old generations
        BA = 0
        do ij = 2, JJ
!          v(-ij+2) = (1d0-damp)*v(-ij+2)+damp*get_W(ij,1)*((util(ij,0)/util(ij,1))**(1d0/egam)-1d0)
            v(-ij+2) = get_W(ij,1)*((util(ij,0)/util(ij,1))**(1d0/egam)-1d0)
            BA(2) = BA(2)+v(-ij+2)*m(ij,1)
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
!          v(it) = (1d0-damp)*v(it) + damp*get_W(1,it)*((ustar/util(1, it))**(1d0/egam)-1d0)
            v(it) = get_W(1,it)*((ustar/util(1, it))**(1d0/egam)-1d0)
            if(it == 2) BA(2) = (BA(2) + v(1))/(1d0+n_p(2))
            if(it > 2) BA(it) = ((1d0+r(it-1))*BA(it-1) + v(it-1))/(1d0+n_p(it))
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, it1, ij

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
        BQ(it) = 0d0
        do ij = 1, JJ
            CC(it) = CC(it)+c(ij,it)*m(ij,it)
            AA(it) = AA(it)+a(ij,it)*m(ij,it)/psi(ij,it)
            LL(it) = LL(it)+h(ij)*m(ij,it)
            BQ(it) = BQ(it)+(1d0-psi(ij,it))*Rn(it)*a(ij,it)*m(ij,it)/psi(ij,it)
        enddo
        AA(it) = AA(it)+a(JJ+1,it)*m(JJ,itm)/(1+n_p(it))
        BQ(it) = BQ(it)+Rn(it)*a(JJ+1,it)*m(JJ,itm)/(1+n_p(it))

        do ij = 1, JJ
            beq(ij, it) = GAM(ij,it)*BQ(it)
!          beq(ij, it) = damp*GAM(ij,it)*BQ(it)+(1d0-damp)*beq(ij,it)
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
        integer :: it1, ij
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


        pen(JR:JJ,it) = kappa(it)*w(it)*h(JR-1)
        ! get budget balancing social security contribution
        Tpen(it) = 0d0
        do ij = JR, JJ
            Tpen(it) = Tpen(it) + pen(ij,it)*m(ij,it)
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
        write(fname,'(a)')'      A      K     BB     BA     BF     BQ       r'
        write(fname,'(6f7.2,f8.2)')AA(it), KK(it), BB(it), BA(it), BF(it), BQ(it), r(it)
        write(fname,'(6f7.2/)')AA(it)/YY(it),KK(it)/YY(it),BB(it)/YY(it),BA(it)/YY(it),BF(it)/YY(it),BQ(it)/YY(it)
        write(fname,'(a)')' Labor Market  '
        write(fname,'(a)')'     LL      w   util     ET'
        write(fname,'(4f7.2/)')LL(it),w(it),util(1,it),ET(it)
        write(fname,'(a)')' GOVERMENT  '
        write(fname,'(a)')'   tauc   tauw   taur   taup   tauk    TXR     DD     rB  '
        write(fname,'(8f7.2)')tauc(it),tauw(it),taur(it),taup(it),tauk(it),TXR(it),&
                           ((1+n_p(it1))*BB(it1)-BB(it)),r(it)*BB(it)
        write(fname,'(35x,3f7.2/)')TXR(it)/YY(it),((1+n_p(it1))*BB(it1)-BB(it))/YY(it),r(it)*BB(it)/YY(it)


        if(fname > 1)then
            write(fname,'(a)')' Age    cons       wn      beq      pen     asset    Diff      '
            do ij = 1, JJ
                it1 = year(it, 1, 2)
                diff = Rn(it)*a(ij,it)+wn(it)*h(ij)+beq(ij,it)+pen(ij,it)-a(ij+1,it1)-p(it)*c(ij,it)
                write(fname,'(i3,6f9.2)')ij,c(ij,it),wn(it)*h(ij),beq(ij,it),pen(ij,it),a(ij,it),diff
            enddo
            write(fname,'(i3,36x,f9.2)')ij,a(JJ+1,it)
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
        write(fname,'(a,a,a)')'         C      A      K     BQ      r', &
            '      w      Y      B     BA     ET   tauc   tauw   taur', &
            '   taup   tauk      v    HEV   y-d'
        diff = YY(0)-CC(0)-GG(0)-II(0)
        ! current generations
        do ij = JJ, 3, -1
            write(fname,'(i3,105x,2f7.2)')-ij+2, v(-ij+2), &
                    ((util(ij,1)/util(ij,0))**(1d0/egam)-1d0)*100d0
        enddo
        write(fname,'(i3,18f7.2)')0,CC(0),AA(0),KK(0),BQ(0),r(0),w(0), &
            YY(0),BB(0),BA(0),ET(0),tauc(0),tauw(0),taur(0),taup(0),tauk(0),v(0), &
            ((util(2,1)/util(2,0))**(1d0/egam)-1d0)*100d0,diff

        ! future generations
        do it = 1, TT
            diff = YY(it)-CC(it)-GG(it)-II(it)
            write(fname,'(i3,18f7.2)')it,CC(it),AA(it),KK(it),BQ(it),r(it),w(it), &
                    YY(it),BB(it),BA(it),ET(it),tauc(it),tauw(it),taur(it),taup(it),tauk(it),v(it), &
                    ((util(1,it)/util(1,0))**(1d0/egam)-1d0)*100d0,diff
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer ::  ik, ik1, iu, ij
        real*8 :: PPs

        util(1,it) = 0d0
        PPs = 1d0
        do ij = 1, JJ
            ik = year(it, 1, ij)
            ik1 = year(it, 1, ij+1)
            PPs = PPs*psi(ij,ik)
            util(1,it) = util(1,it)+beta**(ij-1)*PPs*(c(ij,ik)**egam+&
                        beta*nubeq*(1d0-psi(ij+1,ik1))*(Rn(ik1)*a(ij+1,ik1))**egam)/egam
        enddo

        if(it < 2)then
            do iu = 2, JJ
                util(iu,it) = 0d0
                PPs = 1.0d0
                do ij = iu, JJ
                    ik = year(it, iu, ij)
                    ik1 = year(it, iu, ij+1)
                    util(iu,it) = util(iu,it) + beta**(ij-iu)*PPs*(c(ij,ik)**egam+&
                                beta*nubeq*(1d0-psi(ij+1,ik1))*(Rn(ik1)*a(ij+1,ik1))**egam)/egam
                    if(ij < JJ)PPs = PPs*psi(ij+1,ik1)
                enddo
            enddo
        endif

    end subroutine


    function ET(it)

            implicit none
            integer, intent(in) :: it
            real*8 :: ET, PPs
            integer ::  ik, ik1, ij

            ET = 0d0
            PPs = 1d0
            do ij = 1, JJ
            ik = year(it, 1, ij)
            ik1 = year(it, 1, ij+1)
            PPs = PPs*psi(ij,ik)
            if(ij < JJ)then
                ET = ET + dble(ij)*PPs*(1d0-psi(ij+1,ik1))
            else
                ET = ET + dble(ij)*PPs
            endif
        enddo

    end function

end program
