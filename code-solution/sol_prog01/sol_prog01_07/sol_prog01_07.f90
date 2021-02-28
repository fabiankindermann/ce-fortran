!##############################################################################
! PROGRAM RollDice
!
! ## Simulate the result of rolling n dice with each k sides
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
program RollDice

    use toolbox

    implicit none

    ! declaration of variables and parameters
    integer, parameter :: n = 2
    integer, parameter :: k = 6
    integer, parameter :: iter = 500
    integer ::  i, j, dice(n)
    real*8 :: x, Dsum(n:n*k), SimProb(n:n*k)

    ! set random seed
    call init_random_seed()

    ! simulate dice rolls
    dice = 0
    Dsum = 0d0
    do i = 1, iter
        do j = 1, n
            call random_number(x)
            dice(j) = int(k*x) + 1
        enddo
        Dsum(sum(dice)) = Dsum(sum(dice)) + 1d0
    enddo
    SimProb = Dsum/iter*100d0

    ! print output
    write(*,'(a)')'  Sum  Simulated Probability (in %)'
    do i = n, n*k
        write (*,'(i5,f30.4)')i, SimProb(i)
    enddo

end program
