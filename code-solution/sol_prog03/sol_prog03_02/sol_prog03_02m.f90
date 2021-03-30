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
    real*8, parameter :: Kbar = 10d0
    real*8, parameter :: Lbar = 20d0
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: nu = 0.5d0
    real*8, parameter :: mu = 1d0-1d0/nu
    real*8, parameter :: beta(2) = (/ 0.3d0, 0.6d0/)

    real*8 :: L(2), K(2), Y(2), PP, Ybar, w, r, p(2), U

contains

    ! function to determine market equilibrium
    function markets(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: markets(size(x, 1))

        ! copy prices
        p(1) = 1d0
        p(2) = x(1)
        w    = x(2)
        r    = x(3)

        ! calculate total income
        Ybar = w*Lbar+r*Kbar
        PP   = alpha*p(1)**(1d0-nu)+(1d0-alpha)*p(2)**(1d0-nu)

        ! get market equations
        markets(1) = 1d0/p(1)-(beta(1)/w)**beta(1)*((1d0-beta(1))/r)**(1d0-beta(1))
        markets(2) = 1d0/p(2)-(beta(2)/w)**beta(2)*((1d0-beta(2))/r)**(1d0-beta(2))
        markets(3) = beta(1)*alpha*Ybar/(w*p(1)**(nu-1d0)*PP)+beta(2)*(1-alpha)*Ybar/(w*p(2)**(nu-1d0)*PP)-Lbar

    end function

end module
