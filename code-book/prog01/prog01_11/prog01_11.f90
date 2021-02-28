!##############################################################################
! PROGRAM AlternativeArrays
!
! ## Alternative forms of array value assignments
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
program AlternativeArrays

    implicit none
    integer, parameter :: n = 8
    real*8 :: x(0:n), y(0:n), z(0:n, 0:n)
    integer :: j, k

    ! initialize x
    do j = 0, n
        x(j) = 1d0/dble(n)*dble(j)
    enddo

    ! give y the values of x plus 1
    y(:) = x(:) + 1d0

    ! calculate z
    do j = 0, n
        do k = 0, n
            z(j, k) = x(j)**2 + y(k)
        enddo
    enddo

    ! output table of values
    write(*,'(a,9f7.2)')'    X/Y | ', y(:)
    write(*,'(a)')'        | '
    do j = 0, n
        write(*,'(f7.2,a,9f7.2)')x(j),' | ', z(j, :)
    enddo

end program
