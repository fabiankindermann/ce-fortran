!##############################################################################
! PROGRAM Transport
!
! ## Use simplex algorithm to minimize the total cost of transportation
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program Transport

    use toolbox

    implicit none

    ! number of gravel-pits
    integer, parameter :: n = 3

    ! number of building sites
    integer, parameter :: m = 4

    ! declaration of variables
    real*8 :: c(n*m), x(n*m)
    real*8 :: A(n+m, n*m), b(n+m)

    ! set up transportation costs
    c(:) = (/10d0, 70d0, 100d0, 80d0, & ! gravel-pit 1
            130d0, 90d0, 120d0, 110d0, & ! gravel-pit 2
            50d0, 30d0, 80d0, 10d0/)  ! gravel-pit 3

    ! set up matrix for restriction
    A(1, :) = (/1d0, 1d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/) ! gravel-pit 1
    A(2, :) = (/0d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0/) ! gravel-pit 2
    A(3, :) = (/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 1d0, 1d0, 1d0/) ! gravel-pit 3
    A(4, :) = (/1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0/) ! building site 1
    A(5, :) = (/0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0/) ! building site 2
    A(6, :) = (/0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0/) ! building site 3
    A(7, :) = (/0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0/) ! building site 4

    ! set up target vector for restriction
    b(:) = (/11d0, 13d0, 10d0, 5d0, 7d0, 13d0, 6d0/)

    ! solve the linear problem
    call solve_lin(x, c, A, b, 3, 0, 4)

    ! print output
    write(*,'(a)')'building site  :    1    2    3    4       Sum    costs'
    write(*,'(a,4f5.1,2f10.1)')'gravel-pit 1   :', x(1:4), sum(x(1:4)), sum(x(1:4)*c(1:4))
    write(*,'(a,4f5.1,2f10.1)')'gravel-pit 2   :', x(5:8), sum(x(5:8)), sum(x(5:8)*c(5:8))
    write(*,'(a,4f5.1,2f10.1)')'gravel-pit 1   :', x(9:12), sum(x(9:12)), sum(x(9:12)*c(9:12))
    write(*,'(a,4f5.1,2f10.1)')'Sum            :', x(1) + x(5) + x(9), x(2) + x(6) + x(10), &
        x(3) + x(7) + x(11), x(4) + x(8) + x(12), sum(x), sum(x*c)

end program
