!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
! #VC# VERSION: 1.0  (16 December 2019)
!
!##############################################################################
module globals

    implicit none
    real*8, parameter :: Kbar = 10d0
    real*8, parameter :: Tbar = 30d0
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: alphax = 0.5d0
    real*8, parameter :: nu = 0.5d0
    real*8, parameter :: nux = 0.5d0
    real*8, parameter :: mu = 1d0-1d0/nu
    real*8, parameter :: mux = 1d0-1d0/nux
    real*8, parameter :: beta(2) = (/ 0.3d0, 0.6d0/)

    real*8 :: Ybar, Omega, PP, w, r, p(2), L(2), K(2), Y(2), U
    real*8 :: YD, C, ell

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
        Ybar  = w*Tbar+r*Kbar
        PP    = alphax*p(1)**(1-nux) + (1d0-alphax)*p(2)**(1d0-nux)
        Omega = (1-alpha)*PP**(1-nu) + alpha*w**(1-nu)

        ! get market equations
        markets(1) = 1d0/p(1)-(beta(1)/w)**beta(1)*((1d0-beta(1))/r)**(1d0-beta(1))
        markets(2) = 1d0/p(2)-(beta(2)/w)**beta(2)*((1d0-beta(2))/r)**(1d0-beta(2))
        markets(3) = (beta(1)*alphax/p(1)**(nux-1d0)+beta(2)*(1d0-alphax)/p(2)**(nux-1d0))* &
                     (w**(nu-1d0)*Omega-alpha)/(w**nu*Omega*PP)*Ybar+alpha*Ybar/w**nu/Omega-Tbar

    end function

end module
