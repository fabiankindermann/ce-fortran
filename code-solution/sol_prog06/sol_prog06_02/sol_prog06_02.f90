!##############################################################################
! PROGRAM OLG_SMOPEC
!
! ## The OLG model with transitional dynamics and smopec
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog06_02m.f90"

program OLG_SMOPEC

    use globals
    use toolbox

    implicit none
    integer :: it

    ! initialize variables and government parameters
    call initialize

    ! compute initial long-run equilibrium
    call get_SteadyState

    ! write output
    open(20, file='output.out')
    call output(0, 20)

    ! initialize transitional values
    call get_Transition

    ! write output
    do it = 1, TT
        call output(it, 20)
    enddo
    close(20)

    open(21, file='summary.out')
    call output_summary(21)

    ! get lsra run
    lsra_on = .true.

    ! solve for the transition path
    call get_Transition

    ! write output
    write(21,*)
    call output_summary(21)
    close(21)

end program
