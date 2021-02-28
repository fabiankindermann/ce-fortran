!##############################################################################
! PROGRAM fixedpoint
!
! ## A fixed point iteration scheme to find the root of a function
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
program fixedpoint

    implicit none
    integer :: iter
    real*8 :: xold, x, f, sigma

    ! set initial guess and chose sigma
    xold = 0.05d0
    sigma = 0.2d0

    ! start iteration process
    do iter = 1, 200

        ! calculate function value
        f = 0.5d0*xold**(-0.2d0)+0.5d0*xold**(-0.5d0)-2d0

        ! calculate new value
        x = xold + sigma*f

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
