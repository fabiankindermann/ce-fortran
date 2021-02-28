!##############################################################################
! PROGRAM newton
!
! ## The Newton method to find the root of a function
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
program newton

    implicit none
    integer :: iter
    real*8 :: xold, x, f, fprime

    ! set initial guess
    xold = 0.05d0

    ! start iteration process
    do iter = 1, 200

        ! calculate function value
        f = 0.5d0*xold**(-0.2d0)+0.5d0*xold**(-0.5d0)-2d0

        ! calculate derivative
        fprime = -0.1d0*xold**(-1.2d0)-0.25d0*xold**(-1.5d0)

        ! calculate new value
        x = xold - f/fprime

        write(*,'(i4,f12.7)')iter, abs(x-xold)

        ! check for convergence
        if(abs(x-xold) < 1d-6)then
            write(*,'(/a,f12.7,a,f12.9)')' x = ',x,'    f = ',f
            stop
        endif

        ! copy old value
        xold = x
    enddo

    write(*,'(a)')'Error: no convergence'

end program
