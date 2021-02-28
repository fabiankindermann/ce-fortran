!##############################################################################
! PROGRAM HetAg_Endo
!
! ## The heterogeneous agent model with endogenous gridpoints and distribution
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
include "sol_prog09_13m.f90"

program HetAg_Endo

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
        real*8 :: x
        logical :: check

        ! start timer
        call tic()

        write(*, '(a)')'       K/Y         r      diff_K'

        ! set tolerance level for outer optimization
        call settol_root(sig_out)

        ! set initial guess
        x = r
        check = .false.

        ! solve asset market equilibrium
        call fzero(x, asset_market, check)

        if(check)write(*, '(a)')'ERROR IN ROOTFINDING'

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

        ! initialize grid
        call grid_Cons_Grow(a, a_l, a_u, a_grow)

        ! initialize policy function
        do is = 1, ns
            c(:, is) = 0.04d0*a(:) + eta(is)
            c_new(:, is) = c(:, is)
        enddo

        ! get an initial guess for the distribution function
        phi = 1d0/dble((NA+1)*NS)

    end subroutine


    ! For creating output plots
    subroutine output()

        implicit none

        integer :: ia, is, ial, iar
        real*8 :: CV_c, CV_a, CV_eta, varphi
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, err_temp

        ! aggregate consumption and investment
        CC = sum(phi*c)
        II = delta*KK

        ! calculate coefficients of variation
          CV_eta = 0d0
          CV_c   = 0d0
          CV_a   = 0d0
          do ia = 0, NA
              do is = 1, NS
                  CV_eta = CV_eta + phi(ia, is)*(eta(is)   - LL)**2
                  CV_c   = CV_c   + phi(ia, is)*(c(ia, is) - CC)**2
                  CV_a   = CV_a   + phi(ia, is)*(a(ia)     - KK)**2
              enddo
          enddo
          CV_eta = sqrt(CV_eta)/LL*100d0
          CV_c   = sqrt(CV_c)  /CC*100d0
          CV_a   = sqrt(CV_a)  /AA*100d0

        ! write output to screen
        write(*, '(//a)')'AIYAGARI ECONOMY OUTPUT'
        write(*, '(/a)')'PROD :       K/Y       L/Y         r         w'
        write(*, '(6x,4f10.4)')(/KK, LL/)/YY*100d0, r*100d0, w
        write(*, '(/a)')'GOODS:       C/Y       I/Y       D/Y'
        write(*, '(6x,3f10.4)')(/CC, II, CC+II/)/YY*100d0
        write(*, '(/a)')'               c         a       eta'
        write(*, '(a,4f10.4)')'LEV  :', CC, AA, LL
        write(*, '(a,3f10.4)')'CV   :', CV_c, CV_a, CV_eta
        write(*, '(//)')

        ! plot the distribution of wealth
        call plot(a, sum(phi, 2))
        call execplot(xlim=(/a_l, 30d0/), title='Distribution of Wealth')

        ! plot savings decisions
        do is = 1, NS
            call plot(a, (1d0+r)*a + w*eta(is) - c(:, is))
        enddo
        call execplot(xlim=(/a_l, 0.5d0/), title='Policy Functions')

        ! calculate euler equation error
        call interpolate()
        err = 0d0
        do ia = 0, n_err
            a_com = a_l + (a_u-a_l)*dble(ia)/dble(n_err)
            call linint_Grow(a_com, a_l, a_u, a_grow, NA, ial, iar, varphi)
            do is = 1, NS
                c_err = varphi*c(ial, is) + (1d0-varphi)*c(iar, is)
                is_com = is
                err_temp = abs(foc(c_err)/c_err)
                if(err_temp > err .and. ((1d0+r)*a_com + &
                                  w*eta(is_com) - c_err) > 1d-6)then
                    err = err_temp
                endif
            enddo
        enddo
        write(*, '(a,f15.7)')'Euler equation error:', err

    end subroutine

end program
