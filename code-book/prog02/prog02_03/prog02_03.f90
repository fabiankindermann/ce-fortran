!##############################################################################
! PROGRAM inversion
!
! ## Invert a matrix using LU-decomposition
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
program inversion

    use toolbox

    implicit none
    integer :: i, j
    real*8 :: A(3, 3), Ainv(3, 3), b(3)

    ! set up matrix and vector
    A(1, :) = (/ 2d0,  0d0,  1d0/)
    A(2, :) = (/ 0d0,  4d0,  1d0/)
    A(3, :) = (/ 1d0,  1d0,  2d0/)
    b       = (/30d0, 40d0, 30d0/)

    ! invert A
    Ainv = lu_invert(A)

    ! calculate solution
    b = matmul(Ainv, b)

    ! output
    write(*,'(a,3f7.2/)')'    x = ', (b(j),j=1,3)
    write(*,'(a,3f7.2/,2(8x,3f7.2/))') &
        ' A^-1 = ',((Ainv(i,j),j=1,3),i=1,3)

end program
