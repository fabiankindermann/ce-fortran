!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    implicit none

    ! number of gridpoints to set up intervals
    integer, parameter :: n = 4

    ! number of gridpoints for plot
    integer, parameter :: n_plot = 100

    ! set lower bound for interval
    real*8, parameter :: x_l = 0d0

    ! set upper bound for interval
    real*8, parameter :: x_u = 5d0

    ! set tolerance level
    real*8, parameter :: tol = 1d-6

    ! variables used to locate the global minimum
    real*8 :: x(0:n), minimum_x(n), fmin(n)
    real*8 :: xplot(0:n_plot), yplot(0:n_plot)

contains

    function minimize(a, b)

        implicit none

        ! declaration of variables
        real*8, intent(in) :: a
        real*8, intent(in) :: b
        integer :: iter
        real*8 :: minimize
        real*8 :: a1, b1, x1, x2, f1, f2

        ! set starting values for the iteration
        a1 = a
        b1 = b

        ! start iteration process
        do iter = 1, 200

            ! calculate x1 and x2 values and function values
            x1 = a1 + (3d0-sqrt(5d0))/2d0*(b1-a1)
            x2 = a1 + (sqrt(5d0)-1d0)/2d0*(b1-a1)
            f1 = x1*cos(x1**2)
            f2 = x2*cos(x2**2)

            ! get new values
            if(f1 < f2)then
                b1 = x2
            else
                a1 = x1
            endif

            if(abs(b1-a1) < tol)exit

        enddo

        ! set return value
        if(f1 < f2)then
            minimize = x1
        else
            minimize = x2
        endif

    end function

end module
