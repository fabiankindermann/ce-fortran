!##############################################################################
! PROGRAM RollDiceEndCondition
!
! ## Simulate how long it takes until the game is ended by condition x or y
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
program RollDiceEndCondition

    use toolbox

    implicit none

    ! declaration of variables and parameters
    integer, parameter :: n = 2
    integer, parameter :: mroll = 250
    integer, parameter :: iter = 5000
    integer ::  i, j, k, dice(n), max_roll, x, y
    real*8 :: count_x, SimProb_x

    ! set random seed
    call init_random_seed()

    ! simulate dice rolls
    x = 4
    y = 10
    count_x = 0d0
    max_roll = 0
    do i = 1, iter
        do j = 1, mroll
            do k = 1, n
                call random_int(dice(k), 1, 6)
            enddo

            if(sum(dice) == x)then
                count_x = count_x + 1d0
                max_roll = max(max_roll, j)
                exit
            elseif(sum(dice) == y)then
                max_roll = max(max_roll, j)
                exit
            endif
        enddo
    enddo

    SimProb_x = count_x/iter*100d0

    ! print output
    write(*,'(a, i4)')'Maximum Number of Rolls:', max_roll
    write(*,'(a, f7.2)')'First condition (in %):', SimProb_x
    write(*,'(a, f7.2)')'Second condition (in %):', 100d0 - SimProb_x

contains

    ! simulate a dice roll
    subroutine random_int(result, intl, inth)

        implicit none

        ! declaration of variables
        integer, intent(in) :: intl, inth
        integer, intent(out):: result
        real*8 :: x

        ! call fortran intrinsic reandom number generator
        call random_number(x)
        result = int(inth*x) + intl

    end subroutine

end program
