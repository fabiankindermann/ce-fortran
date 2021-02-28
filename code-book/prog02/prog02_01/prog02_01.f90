!##############################################################################
! PROGRAM matrices
!
! ## Operations on matrices
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
program matrices

    implicit none
    integer :: i, j
    real*8 :: a(4), b(4)
    real*8 :: x(2, 4), y(4, 2), z(2, 2)

    ! initialize vectors and matrices
    a = (/(dble(5-i), i=1, 4)/)
    b = a+4d0
    x(1, :) = (/1d0, 2d0, 3d0, 4d0/)
    x(2, :) = (/5d0, 6d0, 7d0, 8d0/)
    y = transpose(x)
    z = matmul(x, y)

    ! show results of different functions
    write(*,'(a,4f7.1/)')'         vector a = ',(a(i), i=1,4)
    write(*,'(a,f7.1)')  '           sum(a) = ',sum(a)
    write(*,'(a,f7.1/)') '       product(a) = ',product(a)
    write(*,'(a,f7.1)')  '        maxval(a) = ',maxval(a)
    write(*,'(a,i7)')    '        maxloc(a) = ',maxloc(a)
    write(*,'(a,f7.1)')  '        minval(a) = ',minval(a)
    write(*,'(a,i7/)')   '        minloc(a) = ',minloc(a)
    write(*,'(a,4f7.1)') '    cshift(a, -1) = ',cshift(a, -1)
    write(*,'(a,4f7.1/)')'   eoshift(a, -1) = ',eoshift(a, -1)
    write(*,'(a,l7)')    '       all(a<3d0) = ',all(a < 3d0)
    write(*,'(a,l7)')    '       any(a<3d0) = ',any(a < 3d0)
    write(*,'(a,i7/)')   '     count(a<3d0) = ',count(a < 3d0)
    write(*,'(a,4f7.1/)')'         vector b = ',(b(i), i=1,4)
    write(*,'(a,f7.1/)') ' dot_product(a,b) = ',dot_product(a, b)
    write(*,'(a,4f7.1/,20x,4f7.1/)') &
                         '         matrix x = ',((x(i,j), j=1,4), i=1,2)
    write(*,'(a,2f7.1,3(/20x,2f7.1)/)')&
                         '     transpose(x) = ',((y(i,j), j=1,2), i=1,4)
    write(*,'(a,2f7.1/,20x,2f7.1/)') &
                         '      matmul(x,y) = ',((z(i,j), j=1,2), i=1,2)

end program
