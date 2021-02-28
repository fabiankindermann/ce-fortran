!##############################################################################
! PROGRAM golden
!
! ## Golden search for finding the minimum of a function
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
program golden

    implicit none
    real*8, parameter :: p(2) = (/1d0, 2d0/)
    real*8, parameter :: W = 1d0
    real*8 :: a, b, x1, x2, f1, f2
    integer :: iter

    ! initial interval and function values
    a = 0d0
    b = (W-p(1)*0.01d0)/p(2)

    ! start iteration process
    do iter = 1, 200

        ! calculate x1 and x2 and function values
        x1 = a + (3d0-sqrt(5d0))/2d0*(b-a)
        x2 = a + (sqrt(5d0)-1d0)/2d0*(b-a)
        f1 = -(((W-p(2)*x1)/p(1))**0.4d0 + (1d0+x1)**0.5d0)
        f2 = -(((W-p(2)*x2)/p(1))**0.4d0 + (1d0+x2)**0.5d0)

        write(*,'(i4,f12.7)')iter, abs(b-a)

        ! check for convergence
        if(abs(b-a) < 1d-6)then
            write(*,'(/a,f12.7)')' x_1 = ',(W-p(2)*x1)/p(1)
            write(*,'(a,f12.7)')' x_2 = ',x1
            write(*,'(a,f12.7)')' u   = ',-f1
            stop
        endif

        ! get new values
        if(f1 < f2)then
            b = x2
        else
            a = x1
        endif
    enddo

end program
