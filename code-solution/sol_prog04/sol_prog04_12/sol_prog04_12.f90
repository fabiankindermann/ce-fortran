!##############################################################################
! PROGRAM individual_mortality_risk
!
! ## Elimination of individual mortality risk
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
include "sol_prog04_12m.f90"

program individual_mortality_risk

    use toolbox
    use globals

    implicit none

    ! initialize mortality data
    call init_data

    ! initialize interval borders
    call init_borders

    ! simulate individual mortality
    do in = 1, NN

        ! draw random number
        call simulate_uniform(rand, 0d0, 1d0)

        ! compute individual payout and averaging over all contracts
        ind_LI(in) = get_LI(rand)
        agg_LI(in) = (1d0-xi_2)*sum(ind_LI(1:in))/dble(in)

    enddo

    ! write results to the screen
    write(*,'(/a, f10.3/)')'Average Cost = ', agg_LI(NN)

    ! set up and execute plot
    call plot(dble((/(in, in=1, NN)/)), agg_LI)
    call execplot(xlabel='n_i Clients', ylabel='Expected Cost')

contains

    ! initialize mortality and portfolio data
    subroutine init_data()

        implicit none

        !mortality in the base period 0
        q(1) = 0.00186d0
        q(2) = 0.00222d0
        q(3) = 0.00268d0
        q(4) = 0.00389d0
        q(5) = 0.00632d0
        q(6) = 0.01123d0
        q(7) = 0.01957d0
        q(8) = 0.02987d0
        q(9) = 0.04483d0
        q(10) = 0.06528d0
        q(11) = 0.10007d0
        q(12) = 0.17166d0
        q(13) = 0.28709d0
        q(14) = 0.46296d0
        q(15) = 0.66523d0
        q(16) = 1.0d0

    end subroutine init_data

end program
