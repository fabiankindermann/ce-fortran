!##############################################################################
! PROGRAM HetAg_Lab
!
! ## Heterogeneous Agent model with variable labor supply
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
include "sol_prog09_14m.f90"

program HetAg_Lab

    use globals

    implicit none

    ! set initial guess for interest rate
    r = 0.04d0

    ! initialize variables
    call initialize()

    ! solve initial equilibrium
    call solve_model()

contains


    ! For solving the model
    subroutine solve_model

        use toolbox

        implicit none
        real*8 :: x
        logical :: check

        ! start timer
        call tic()

        write (*,'(a)')'       K/Y         r      diff_K'

        ! set tolerance level for outer optimization
        call settol_root(sig_out)

        ! set initial guess
        x = KK/LL
        check = .false.

        ! solve asset market equilibrium
        call fzero(x, asset_market, check)

        if(check)write (*,'(a)')'ERROR IN ROOTFINDING'

        call toc

        ! write output to the screen
        call output

    end subroutine


    ! For initializing variables
    subroutine initialize()

        use toolbox

        implicit none
        integer :: is, ia
        real*8 :: temp

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi, weights)
        eta = exp(eta)

        ! calculate aggregate labor supply and capital
        LL = sum(eta*weights)
        KK = LL*(alpha/(r+delta))**(1d0/(1d0-alpha))

        ! initialize grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize policy function
        w = (1d0-alpha)*(KK/LL)**alpha
        do is = 1, NS
            do ia = 0, NA
                aplus(ia, is) = a(ia)
                l(ia, is) = min(max(nu+(1d0-nu)*(aplus(ia, is)-(1d0+r)*a(ia))/(w*eta(is)), 0d0), 1d0-1d-10)
                c(ia, is) = max((1d0+r)*a(ia)+w*eta(is)*l(ia, is)-aplus(ia, is), 1d-10)
            enddo
        enddo

        ! get an initial guess for the distribution function
        phi = 1d0/dble((NA+1)*NS)

    end subroutine


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer :: ia, is
        real*8 :: CV_c, CV_a, CV_l

        ! aggregate consumption and investment
        CC = sum(phi*c)
        II = delta*KK

        ! calculate coefficients of variation
        CV_c = 0d0
        CV_a = 0d0
        CV_l = 0d0
        do ia = 0, NA
            do is = 1, NS
                CV_c   = CV_c   + phi(ia, is)*(c(ia, is)    - CC)**2
                CV_a   = CV_a   + phi(ia, is)*(a(ia)        - KK)**2
                CV_l   = CV_l   + phi(ia, is)*(leff(ia, is) - LL)**2
            enddo
        enddo
        CV_c = sqrt(CV_c)/CC*100d0
        CV_a = sqrt(CV_a)/AA*100d0
        CV_l = sqrt(CV_l)/LL*100d0

        ! write output to screen
        write(*,'(//a)')'AIYAGARI ECONOMY OUTPUT'
        write(*,'(/a)')'PROD :       K/Y       L/Y         r         w'
        write(*,'(6x,4f10.4)')(/KK, LL/)/YY*100d0, r*100d0, w
        write(*,'(/a)')'GOODS:       C/Y       I/Y       D/Y'
        write(*,'(6x,3f10.4)')(/CC, II, CC+II/)/YY*100d0
        write(*,'(/a)')'               c         a         l'
        write(*,'(a,4f10.4)')'LEV  :', CC, AA, LL
        write(*,'(a,3f10.4)')'CV   :', CV_c, CV_a, CV_l
        write(*,'(//)')

        ! plot the distribution of wealth
        call plot(a, sum(phi, 2))
        call execplot(xlim=(/a_l, 30d0/), title='Distribution of Wealth')

        ! plot savings decisions
        do is = 1, NS
            call plot(a, (1d0+r)*a + w*eta(is) - c(:, is))
        enddo
        call execplot(xlim=(/a_l, 0.5d0/), title='Savings Decisions')

        do is = 1, NS
            call plot(a, l(:, is))
        enddo
        call execplot(xlim=(/a_l, 0.5d0/), title='Labor Supply')

    end subroutine

end program
