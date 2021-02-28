!##############################################################################
! PROGRAM Fibonacci
!
! ## Compute the n-th element of the Fibonacci-Series
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
program Fibonacci

    use toolbox

    implicit none

    ! declaration of parameters
    integer, parameter :: plotmax = 100
    integer :: n, i
    real*8 :: x(1:plotmax), diff(1:plotmax)

    ! read input
    write(*,'(a)')'Type in the requested element of series: '
    read(*,*)n

    ! print outut
    write(*,'(i3,a)') n, '-element of the Fibonacci series (with do-loops):           '
    write(*,'(f40.12)') fib(n)
    write(*,'(i3,a)') n, '-element of the Fibonacci series (with Binets formula): '
    write(*,'(f40.12)') Binform(n)

    ! plot the difference between iterative solution Binet's-formula
    do i = 1, plotmax
        x(i) = i
        diff(i) = abs(fib(i) - Binform(i))
    enddo

    call plot(x, diff, legend='Diff')
    call execplot()

contains

    ! iterative solution
    function fib(n)

        implicit none

        ! input variables
        integer, intent(in):: n
        integer:: j
        real*8 :: a0, a1, a2, fib

        ! Calculation of the result
        if(n <= 2)then
            a0 = 1
        else

            ! initialize the first two elements
            a1 = 1
            a2 = 1

            ! calculate the n-th element
            do j = 3, n
                a0 = a1 + a2
                a2 = a1
                a1 = a0
            enddo
        endif

        fib = a0

    end function

    ! closed-form solution
    function Binform(n)

        implicit none

        ! input variables
        integer, intent(in):: n
        real*8:: Binform
        real*8:: golden_ratio

        golden_ratio = (1d0 + sqrt(5d0))/2d0

        Binform = (golden_ratio**(n) - (1d0 - golden_ratio)**(n))/sqrt(5d0)

    end function

end program
