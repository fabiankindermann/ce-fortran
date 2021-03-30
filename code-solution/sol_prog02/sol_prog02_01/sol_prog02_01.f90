!##############################################################################
! PROGRAM Matrices
!
! ## Perform lu-decomposition and solve linear equation system
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program Matrices

    use toolbox

    implicit none

    ! declaration of variables
    real*8 :: A(4, 4), L(4, 4), U(4, 4), A_test(4, 4), b(4), b_test(4)
    integer :: i, j

    ! initialize matrix A
    A(1, :) = (/1d0, 5d0, 2d0, 3d0/)
    A(2, :) = (/1d0, 6d0, 8d0, 6d0/)
    A(3, :) = (/1d0, 6d0, 11d0, 2d0/)
    A(4, :) = (/1d0, 7d0, 17d0, 4d0/)

    ! initialize vector b
    b = (/1d0, 2d0, 1d0, 1d0/)

    ! print matrix A to the screen
    write(*,'(a)')'A = '
    write(*,'(4f7.1)')((A(j, i),i = 1, 4), j = 1, 4)

    ! decompose matrix A
    call lu_dec(A, L, U)

    ! check the result from lu_dec
    A_test = matmul(L, U)

    ! print output
    write(*,'(a)')'L = '
    write(*,'(4f7.1)')((L(j, i), i = 1,4), j = 1, 4)
    write(*,'(/a)')'U = '
    write(*,'(4f7.1)')((U(j, i), i = 1, 4), j = 1, 4)
    write(*,'(/a)')'A_test = '
    write(*,'(4f7.1)')((A_test(j, i), i = 1, 4), j = 1, 4)
    write(*,'(/a/)')'-----------------------------'

    ! solve the linear equation system
    call lu_solve(A, b)

    ! check the result from lu_solve
    b_test = matmul(A, b)

    ! print output
    write(*,'(a)')'b = '
    write(*,'(4f7.1)')(b(i), i = 1, 4)
    write(*,'(/a)')'b_test = '
    write(*,'(4f7.1)')(b_test(i), i = 1, 4)

end program
