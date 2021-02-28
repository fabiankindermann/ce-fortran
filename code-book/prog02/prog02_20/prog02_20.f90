!##############################################################################
! PROGRAM simplex_alg
!
! ## The simplex algorithm for solving linear constrained minimization problems
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
program simplex_alg

    use toolbox

    implicit none
    real*8 :: c(2), x(2), A(3, 2), b(3)
    integer :: j

    ! set up matrix, target vector and coefficients
    A(1, :) = (/ 1d0,  1d0/)
    A(2, :) = (/ 4d0,  1d0/)
    A(3, :) = (/20d0, 10d0/)

    b(:) = (/100d0, 160d0, 1100d0/)

    c(:) = (/-120d0, -40d0/)

    ! solve linear program
    call solve_lin(x, c, A, b, 3, 0, 0)

    ! output
    write(*,'(/a,2f10.2)')'      x = ',(x(j), j=1,2)
    write(*,'(/a,f10.2)')' Cons 1 = ', b(1)-sum(A(1, :)*x)
    write(*,'(a,f10.2)')' Cons 2 = ', b(2)-sum(A(2, :)*x)
    write(*,'(a,f10.2)')' Cons 3 = ', b(3)-sum(A(3, :)*x)

end program
