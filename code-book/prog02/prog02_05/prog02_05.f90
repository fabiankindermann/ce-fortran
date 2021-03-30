!##############################################################################
! PROGRAM bisection
!
! ## A bisection method to find the root of a function
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program bisection

    implicit none
    integer :: iter
    real*8 :: x, a, b, fx, fa, fb

    ! set initial guesses and function values
    a =  0.05d0
    b =  0.25d0
    fa = 0.5d0*a**(-0.2d0)+0.5d0*a**(-0.5d0)-2d0
    fb = 0.5d0*b**(-0.2d0)+0.5d0*b**(-0.5d0)-2d0

    ! check whether there is a root in [a,b]
    if(fa*fb >= 0d0)then
        stop 'Error: There is no root in [a,b]'
    endif

    ! start iteration process
    do iter = 1, 200

        ! calculate new bisection point and function value
        x = (a+b)/2d0
        fx = 0.5d0*x**(-0.2d0)+0.5d0*x**(-0.5d0)-2d0

        write(*,'(i4,f12.7)')iter, abs(x-a)

        ! check for convergence
        if(abs(x-a) < 1d-6)then
            write(*,'(/a,f12.7,a,f12.9)')' x = ',x,'    f = ',fx
            stop
        endif

        ! calculate new interval
        if(fa*fx < 0d0)then
            b = x
            fb = fx
        else
            a = x
            fa = fx
        endif
    enddo

    write(*,'(a)')'Error: no convergence'

end program
