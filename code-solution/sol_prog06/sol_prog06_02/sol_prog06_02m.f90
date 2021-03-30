!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    implicit none

    ! model parameters
    integer, parameter :: TT = 25
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.9d0
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: delta = 0d0
    real*8, parameter :: tol = 1d-5
    real*8, parameter :: damp = 0.25d0
    integer, parameter :: itermax = 1000

    ! model variables
    real*8 :: g(3), by(0:TT), kappa(0:TT), n_p(0:TT)
    integer :: tax(0:TT)
    logical :: lsra_on, smopec
    real*8 :: w(0:TT), r(0:TT), wn(0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT), pen(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BA(0:TT), BF(0:TT), TB(0:TT)
    real*8 :: a(3,0:TT), c(3,0:TT), util(3,0:TT), v(-1:TT)

contains


    ! initialize variables and government parameters
    subroutine initialize

        ! set baseline parameters
        g(1:3)      = (/0.12d0, 0.12d0, 0.0d0/)
        by(0:TT)    = 0d0
        kappa(0:TT) = 0d0
        n_p(0:TT)   = 0.2d0
        tax(0:TT)   = 1
        lsra_on     = .false.
        smopec      = .true.

        ! set reform values (uncomment respective line for different tables)
        tax(1:TT) = 2                               ! Table 6.2
        !tax(1:TT) = 4                               ! Table 6.3
        !tax(1:TT) = 3                               ! Table 6.4
        !tax(1:TT) = 3 ; by(1:TT) = -0.058857d0      ! Table 6.5
        !kappa(1:TT) = 0.5d0                         ! Table 6.6
        !by(1:TT) = 0.0986d0                         ! Table 6.7
        !n_p(1:TT) = 0.0d0                           ! Table 6.8
        !n_p(1:TT) = 0.0d0 ; kappa = 0.5d0 ; g = 0d0 ! Table 6.9

        ! get labor supply and pension payments
        LL = (2d0+n_p)/(1d0+n_p)
        taup = kappa/((2d0+n_p)*(1d0+n_p))

        ! initialize tax rates
        tauc = 0d0
        tauw = 0d0
        taur = 0d0

        ! initialize assets, LSRA payments, debt holdings and foreign variables
        a = 0d0
        v = 0d0
        BA = 0d0
        TB = 0d0
        BF = 0d0

    end subroutine


    ! solves for steady state using Gauss-Seidel
    subroutine get_SteadyState

        implicit none
        integer :: iter

        write(*,'(a/)')'INITIAL EQUILIBRIUM'

        ! initial guess for capital stock
        KK(0) = 1d0

        do iter = 1, itermax

            ! get prices, decisions and quantities
            call factor_prices(0)
            call decisions(0)
            call quantities(0)
            call government(0)

            ! print output to screen
            write(*,'(a,i4,a,f16.10)')'Iteration: ',iter,' Diff: ', &
                abs((YY(0) - CC(0) - II(0) - GG(0) - TB(0))/YY(0))

            ! check whether goods market is in equilibrium
            if(abs(YY(0) - CC(0) - II(0) - GG(0) - TB(0))/YY(0) < tol) &
                    exit

        enddo

        if(iter > itermax)then
            write(*, '(/a/)') '!!! No equilibrium found !!!'
        else
            write(*,*)
        endif

    end subroutine


    ! solves for transition path using Gauss-Seidel
    subroutine get_Transition

        implicit none
        integer :: it, nmarket, iter

        if(lsra_on)then
            write(*,'(/a/)')'TRANSITION PATH WITH LSRA'
        else
            write(*,'(/a/)')'TRANSITION PATH WITHOUT LSRA'
        endif

        ! initialize values from initial equilibrium
        KK(:) = KK(0)
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
                if(abs(YY(it) - CC(it) - II(it) - GG(it) - TB(it))/YY(it) < tol) &
                    nmarket = nmarket + 1
            enddo

            write(*,'(a,i4,a,i4,a,f16.10)')'Iteration: ',iter, &
                ' Markets: ', nmarket,' Diff: ', &
                maxval(abs(YY(:) - CC(:) - II(:) - GG(:) - TB(:))/YY(:))
            if(nmarket == TT)exit

        enddo

        if(iter > itermax)then
            write(*,'(/a/)') '!!! No equilibrium found !!!'
        else
            write(*,*)
        endif

    end subroutine


    ! subroutine for calculating factor prices in a certain year
    subroutine factor_prices(it)

        implicit none
        integer, intent(in) :: it

        ! interest payments and wages in smopec for t > 0
        if(smopec .and. it > 0)then
            r(it) = r(0)
        else
            r(it) = alpha*(KK(it)/LL(it))**(alpha-1d0)-delta
        endif

        ! factor prices and pension payments in year t
        w(it) = (1d0-alpha)*(KK(it)/LL(it))**alpha
        wn(it) = w(it)*(1d0-tauw(it)-taup(it))
        Rn(it) = 1d0+r(it)*(1d0-taur(it))
        p(it) = 1d0 + tauc(it)
        pen(it) = kappa(it)*w(max(it-1, 0))

    end subroutine


    ! subroutine for calculating individual decisions in a certain year
    subroutine decisions(it)

        implicit none
        integer, intent(in) :: it
        integer :: it1, it2, itm
        real*8 :: PVI, PSI

        ! calculate future and past years
        it1 = year(it, 1, 2)
        it2 = year(it, 1, 3)
        itm = year(it, 2, 1)

        ! individual decisions
        PVI = wn(it) + wn(it1)/Rn(it1) + pen(it2)/(Rn(it1)*Rn(it2)) + v(it)
        PSI = p(it)*(1d0 + beta**gamma*(p(it1)/p(it)/Rn(it1))**(1d0-gamma) + &
              beta**(2*gamma)*(p(it2)/p(it)/Rn(it1)/Rn(it2))**(1d0-gamma))
        c(1, it) = PVI/PSI

        if(it == 1)then
            PVI = Rn(it)*a(2, 0) + wn(it) + pen(it1)/Rn(it1) + v(0)
            PSI = p(it)*(1d0 + beta**gamma*(p(it1)/p(it)/Rn(it1))**(1d0-gamma))
            c(2, it) = PVI/PSI
            c(3, it) = (pen(it) + Rn(it)*a(3, itm) + v(-1))/p(it)
            a(2, it) = wn(itm) - p(itm)*c(1, itm)
        else
            c(2, it) = (beta*Rn(it)*p(itm)/p(it))**gamma*c(1, itm)
            c(3, it) = (beta*Rn(it)*p(itm)/p(it))**gamma*c(2, itm)
            a(2, it) = wn(itm) + v(itm) - p(itm)*c(1, itm)
        endif

        if(it == 2)then
            a(3, it) = wn(itm) + Rn(itm)*a(2, itm) + v(0) - p(itm)*c(2, itm)
        else
            a(3, it) = wn(itm) + Rn(itm)*a(2, itm) - p(itm)*c(2, itm)
        endif

    end subroutine


    ! calculates year at which age ij agent is ij_p
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
        integer :: it, it1, it2
        real*8 :: PVI(-1:TT), sum1, sum2, ustar, PVV

        ! calculate utility for each generation
        do it = 1, TT
            call utility(it)
        enddo

        ! transfers to old generations
        PVI(-1) = Rn(1)*a(3,1) + pen(1) + v(-1)
        PVI(0)  = Rn(1)*a(2,1) + wn(1) + pen(2)/Rn(2) + v(0)
        v(-1) = v(-1) + PVI(-1)*((util(3,0)/util(3,1))**(1d0/egam)-1d0)
        v(0)  = v(0) + PVI(0) *((util(2,0)/util(2,1))**(1d0/egam)-1d0)
        BA(2) = v(-1)/((1d0+n_p(0))*(1d0+n_p(1))) + v(0)/(1d0+n_p(1))

        ! long run equilibrium
        PVI(TT) = wn(TT) + wn(TT)/Rn(TT) + pen(TT)/Rn(TT)**2 + v(TT)
        PVV  = v(TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum1 = PVI(TT)*(1d0+r(TT))/(r(TT)-n_p(TT))
        sum2 = PVI(TT)*(util(1,TT)*egam)**(-1d0/egam)*(1d0+r(TT))/ &
               (r(TT)-n_p(TT))

        ! transition path
        do it = TT-1, 1, -1
            it1 = year(it, 1, 2)
            it2 = year(it, 1, 3)
            PVI(it) = wn(it) + wn(it1)/Rn(it1) + pen(it2)/(Rn(it1)*Rn(it2)) + v(it)
            PVV = PVV*(1d0+n_p(it1))/(1d0+r(it1)) + v(it)
            sum1 = sum1*(1d0+n_p(it1))/(1d0+r(it1)) + PVI(it)
            sum2 = sum2*(1d0+n_p(it1))/(1d0+r(it1)) + PVI(it)* &
                   (util(1,it)*egam)**(-1d0/egam)
        enddo

        ! calculate ustar for future generations
        ustar = ((sum1-BA(2)-PVV)/sum2)**egam/egam

        ! calculate transfers to future generations and debt of LSRA
        do it = 1, TT
            v(it) = v(it) + PVI(it)*((ustar/util(1, it))**(1d0/egam)-1d0)
            if(it == 2)BA(2) = (BA(2) + v(1))/(1d0+n_p(2))
            if(it > 2)BA(it) = ((1d0+r(it-1))*BA(it-1) + v(it-1))/(1d0+n_p(it))
        enddo

    end subroutine


    ! subroutine for calculating quantities in a certain year
    subroutine quantities(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm, itp

        itm = year(it, 2, 1)
        itp = year(it, 1, 2)

        ! individual decisions
        YY(it) = KK(it)**alpha * LL(it)**(1d0-alpha)
        CC(it) = c(1, it) + c(2, it)/(1d0+n_p(it)) + &
            c(3, it)/((1d0+n_p(it))*(1d0+n_p(itm)))
        GG(it) = g(1) + g(2)/(1d0+n_p(it)) + &
                 g(3)/((1d0+n_p(it))*(1d0+n_p(itm)))
        AA(it) = a(2, it)/(1d0+n_p(it)) + a(3, it)/((1d0+n_p(it))*(1d0+n_p(itm)))
        BB(it) = by(itm)*YY(it)

        ! small open or closed economy
        if(smopec .and. it > 0)then
            KK(it) = LL(it)*((r(it)+delta)/alpha)**(1d0/(alpha-1d0))
            BF(it) = AA(it) - KK(it) - BA(it) - BB(it)
            TB(it) = (1d0+n_p(itp))*BF(itp) - (1d0+r(it))*BF(it)
        else
            KK(it) = damp*(AA(it)-BB(it)-BA(it)) + (1d0-damp)*KK(it)
        endif

        II(it) = (1d0+n_p(itp))*KK(itp)-(1d0-delta)*KK(it)

    end subroutine


    ! subroutine for calculating government parameters
    subroutine government(it)

        implicit none
        integer, intent(in) :: it
        integer :: itm

        itm = year(it, 2, 1)

        ! get budget balancing tax rate
        if(tax(it) == 1)then
            tauc(it) = ((1d0+r(it))*BB(it) + GG(it) - (tauw(it)*w(it)*LL(it) + &
                taur(it)*r(it)*AA(it) + (1d0+n_p(it))*BB(min(it+1,TT))))/CC(it)
        elseif(tax(it) == 2)then
            tauw(it) = ((1d0+r(it))*BB(it) + GG(it) - (tauc(it)*CC(it) + &
                (1d0+n_p(it))*BB(min(it+1,TT))))/(w(it)*LL(it) + r(it)*AA(it))
            taur(it) = tauw(it)
        elseif(tax(it) == 3)then
            tauw(it) = ((1d0+r(it))*BB(it) + GG(it) - (tauc(it)*CC(it) + &
                taur(it)*r(it)*AA(it) + (1d0+n_p(it))*BB(min(it+1,TT))))/(w(it)*LL(it))
        else
            taur(it) = ((1d0+r(it))*BB(it) + GG(it) - (tauc(it)*CC(it) + &
                tauw(it)*w(it)*LL(it) + (1d0+n_p(it))*BB(min(it+1,TT))))/(r(it)*AA(it))
        endif

        ! get budget balancing social security contribution
        taup(it) = (pen(it)/((2d0+n_p(it))*(1d0+n_p(itm))))/w(it)

    end subroutine


    ! subroutine for writing output to a file
    subroutine output(it, fname)

        implicit none
        integer, intent(in) :: it
        integer, intent(in) :: fname
        real*8 :: diff

        ! household utility
        call utility(it)
        II(0) = (n_p(0) + delta)*KK(0)

        ! difference on goods market
        diff = YY(it) - CC(it) - GG(it) - II(it) - TB(it)

        ! output
        write(fname,'(a,i3/)')'  Equilibrium: Year', it
        write(fname,'(a)')'     c1     c2     c3      Y      w      r      U '
        write(fname,'(7f7.2/)')c(:, it), YY(it), w(it), r(it), util(1, it)

        write(fname,'(a)')'     a2     a3      K      L     BF'
        write(fname,'(5f7.2/)')a(2, it), a(3, it), KK(it), LL(it), BF(it)

        write(fname,'(a)')'   tauw   taur   tauc   taup    pen      B'
        write(fname,'(6f7.2/)')tauw(it), taur(it), tauc(it), taup(it), pen(it), BB(it)

        write(fname,'(a)')'      Y      C      G      I     TB      DIFF'
        write(fname,'(5f7.2,f10.4//)')YY(it), CC(it), GG(it), II(it), TB(it), diff

    end subroutine


    ! subroutine for writing output summary to a file
    subroutine output_summary(fname)

        implicit none
        integer, intent(in) :: fname
        integer :: it
        real*8 :: diff

        ! headline
        write(fname,'(a,a,a)')'        c1     c2     c3      A      K      r', &
            '      w      Y      B     BA     BF     TB   tauc   tauw   taur', &
            '   taup      v    HEV   y-d'

        ! current generations
        write(fname,'(i3,112x,2f7.2)')-1, v(-1), &
            ((util(3, 1)/util(3, 0))**(1d0/egam)-1d0)*100d0
        write(fname,'(i3,18f7.2)')0, c(1, 0), c(2, 0), c(3, 0), AA(0), KK(0), r(0), w(0), &
            YY(0), BB(0), BA(0), BF(0), TB(0), tauc(0), tauw(0), taur(0), taup(0), v(0), &
            ((util(2, 1)/util(2, 0))**(1d0/egam)-1d0)*100d0

        ! future generations
        do it = 1, TT
            diff = YY(it) - CC(it) - GG(it) - II(it) - TB(it)
            write(fname,'(i3,18f7.2,f7.3)')it, c(1, it), c(2, it), c(3, it), AA(it), &
                KK(it), r(it), w(it), YY(it), BB(it), BA(it), BF(it), TB(it), tauc(it), tauw(it), &
                taur(it), taup(it), v(it), &
                ((util(1, it)/util(1, 0))**(1d0/egam)-1d0)*100d0, diff
        enddo

    end subroutine


    ! subroutine to compute household utility
    subroutine utility(it)

        implicit none
        integer, intent(in) :: it
        integer :: it1, it2

        ! get future years
        if(it == 0)then
            it1 = 0
            it2 = 0
        else
            it1 = min(it+1, TT)
            it2 = min(it+2, TT)
        endif

        ! oldest cohort
        util(3, it) = c(3, it)**egam/egam

        ! middle cohort
        util(2, it) = c(2, it)**egam/egam + beta*c(3, it1)**egam/egam

        ! youngest cohort
        util(1, it) = c(1, it)**egam/egam + beta*c(2, it1)**egam/egam + &
                                            beta**2*c(3, it2)**egam/egam


    end subroutine

end module
