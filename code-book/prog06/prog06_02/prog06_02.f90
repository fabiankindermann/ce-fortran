!##############################################################################
! PROGRAM TR_OLG
!
! ## The OLG model with transitional dynamics
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
include "prog06_02m.f90"

program TR_OLG

    use globals
    use toolbox

    implicit none
    logical :: check
    real*8 :: x(2)
    integer :: it

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

    ! initialize assets, LSRA payments and debt holdings
    a = 0d0
    v = 0d0
    BA = 0d0

    ! get initial guess
    x(:) = 0.7d0

    ! solve the steady state equation system
    call fzero(x, eqns_Initial, check)

    ! check whether the solution is valid
    if(check)then
        write(*, '(/a/)') '!!! No equilibirium found !!!'
        stop
    endif

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
