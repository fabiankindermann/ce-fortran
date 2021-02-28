!##############################################################################
! PROGRAM jacobi
!
! ## Use Jacobi iteration algorithm to solve linear equation system
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
program jacobi

    implicit none
    integer :: iter, i
    real*8 :: A(3, 3), Dinv(3, 3), ID(3, 3), C(3, 3)
    real*8 :: b(3), d(3), x(3), xold(3)

    ! set up matrices and vectors
    A(1, :) = (/ 2d0,  0d0,  1d0/)
    A(2, :) = (/ 0d0,  4d0,  1d0/)
    A(3, :) = (/ 1d0,  1d0,  2d0/)
    b       = (/30d0, 40d0, 30d0/)

    ID = 0d0
    Dinv = 0d0
    do i = 1, 3
        ID(i, i) = 1d0
        Dinv(i, i) = 1d0/A(i, i)
    enddo

    ! calculate iteration matrix and vector
    C = ID-matmul(Dinv, A)
    d = matmul(Dinv, b)

    ! initialize xold
    xold = 0d0

    ! start iteration
    do iter = 1, 200

        x = d + matmul(C, xold)

        write(*,'(i4,f12.7)')iter, maxval(abs(x-xold))

        ! check for convergence
        if(maxval(abs(x-xold)) < 1d-6)then
            write(*,'(/a,3f12.2)')' x = ',(x(i),i=1,3)
            stop
        endif

        xold = x
    enddo

    write(*,'(a)')'Error: no convergence'

end program
