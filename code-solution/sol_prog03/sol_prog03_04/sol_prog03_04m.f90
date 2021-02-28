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
    real*8, parameter :: sigma(2) = (/ 0.5d0, 0.5d0/)
    real*8, parameter :: rho(2) = 1d0-1d0/sigma

    real*8 :: L(2), K(2), Y(2), ly(2), ky(2)
    real*8 :: Ybar, Omega, PP, w, r, p(2), U
    real*8 :: YD, C, ell

contains

    ! function to determine market equilibrium
    function markets(x)

        implicit none
        real*8, intent(in) :: x
        real*8 :: markets

        ! calculate prices
        w     = 1d0
        r     = x
        ly    = (beta + (1d0-beta)*(beta*r/((1d0-beta)*w))**(1d0-sigma))**(-1d0/rho)
        ky    = ((1d0-beta) + beta*((1d0-beta)*w/(beta*r))**(1d0-sigma))**(-1d0/rho)
        p     = ly*w + ky*r
        PP    = alphax*p(1)**(1d0-nux) + (1d0-alphax)*p(2)**(1d0-nux)
        Omega = (1d0-alpha)*PP**(1d0-nu) + alpha*w**(1d0-nu)

        ! calculate other economic variables
        Ybar = w*Tbar+r*Kbar
        C    = (1d0-alpha)*Ybar/PP**nu/Omega
        ell  = alpha*Ybar/w**nu/Omega
        YD   = Ybar-w*ell
        Y(1) = alphax*YD/p(1)**nux/PP
        Y(2) = (1d0-alphax)*YD/p(2)**nux/PP
        L    = ly*Y
        K    = ky*Y

        ! get market equations
        markets = K(1)+K(2)-Kbar

    end function

end module
