!##############################################################################
! PROGRAM HetAg_Sim
!
! ## The heterogeneous agent model with simulation methods
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
include "prog09_08m.f90"

program HetAg_Sim

    use globals

    implicit none

    ! initialize variables
    call initialize()

    ! set initial guess of the interest rate
    r = 0.040d0

    ! solve initial equilibrium
    call solve_model()

contains


    ! For solving the model
    subroutine solve_model

        use toolbox

        implicit none
        real*8 :: x_in
        logical :: check

        ! start timer
        call tic()

        write(*, '(a)')'       K/Y         r      diff_K'

        ! set tolerance level for outer optimization
        call settol_root(sig_out)

        ! set initial guess
        x_in = r
        check = .false.

        ! solve asset market equilibrium
        call fzero(x_in, asset_market, check)

        if(check)write(*,'(a)')'ERROR IN ROOTFINDING'

        call toc

        ! write output to the screen
        call output

    end subroutine


    ! For initializing variables
    subroutine initialize()

        use toolbox

        implicit none
        integer :: is

        ! calculate the shock process
        call discretize_AR(rho, 0d0, sigma_eps, eta, pi, weights)
        eta = exp(eta)

        ! calculate aggregate labor supply
        LL = sum(eta*weights)

        ! simulate draws from the AR(1) process
        call simulate_AR(pi, is_t)

        ! initialize grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize policy function
        do is = 1, ns
            c(:, is) = 0.04d0*a(:) + eta(is)
            c_new(:, is) = c(:, is)
        enddo

    end subroutine


    ! For creating output plots
    subroutine output()

        implicit none
        real*8 :: CV_c, CV_a, CV_eta

        ! aggregate consumption and investment
        CC = sum(c_t(0:TT))/dble(TT+1)
        II = delta*KK

        ! calculate coefficients of variation
        CV_eta = sqrt(sum((eta_t(0:TT)-LL)**2)/dble(TT))/LL*100d0
        CV_c   = sqrt(sum((c_t(0:TT)-CC)**2)/dble(TT))/CC*100d0
        CV_a   = sqrt(sum((a_t(0:TT)-KK)**2)/dble(TT))/AA*100d0

        ! write output to screen
        write(*, '(//a)')'AIYAGARI ECONOMY OUTPUT'
        write(*,'(/a)')'PROD :       K/Y       L/Y         r         w'
        write(*,'(6x,4f10.4)')(/KK, LL/)/YY*100d0, r*100d0, w
        write(*,'(/a)')'GOODS:       C/Y       I/Y       D/Y'
        write(*,'(6x,3f10.4)')(/CC, II, CC+II/)/YY*100d0
        write(*,'(/a)')'               c         a       eta'
        write(*,'(a,4f10.4)')'LEV  :', CC, AA, LL
        write(*,'(a,3f10.4)')'CV   :', CV_c, CV_a, CV_eta
        write(*,'(//)')

    end subroutine

end program
