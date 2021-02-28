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
    real*8, parameter :: alpha(2) = (/0.3d0, 0.4d0/)
    real*8, parameter :: beta(2) = (/0.3d0, 0.6d0/)
    real*8, parameter :: a0(2) = (/0.2d0, 0.2d0/)
    real*8, parameter :: G = 3d0
    real*8, parameter :: ID(2, 2) = reshape((/1d0, 0d0, 0d0, 1d0/),(/2,2/))
    real*8            :: a(2, 2) = reshape((/0d0, 0.3d0, 0.2d0, 0d0/),(/2,2/))
    real*8            :: tauw = 0d0
    real*8            :: taur = 0d0
    real*8            :: tau(2) = 0d0

    integer :: lambda
    real*8 :: L(2), K(2), Y(2), ly(2), ky(2), Xd(2), T(2)
    real*8 :: Ybarn, w, wn, r, rn, p(2), q(2), U, ell

contains

    ! function to determine market equilibrium
    function markets(x)

        use toolbox

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: markets(size(x, 1))

        ! copy producer prices and taxes
        w      = 1d0
        r      = x(1)
        tau(1) = x(2)
        tau(2) = 0d0

        ! 1. calculate K/Y and L/Y
        ky = a0*((1d0-beta)/beta*w/r)**beta
        ly = a0*(beta/(1d0-beta)*r/w)**(1d0-beta)

        ! 2. determine producer prices
        q = w*ly+r*ky
        a(1,2) = a(1,2)*(1d0+lambda*tau(1))
        call lu_solve(ID-transpose(a), q)

        ! 3. consumer prices and demands
        p  = q*(1d0+tau)
        wn = w*(1d0-tauw)
        rn = r*(1d0-taur)
        Ybarn = wn*Tbar+rn*Kbar
        Xd = alpha/p*Ybarn
        ell  = (1d0-alpha(1)-alpha(2))/wn*Ybarn

        ! 4. determine output levels
        Y(1) = Xd(1)+G
        Y(2) = Xd(2)
        a(1,2) = a(1,2)/(1d0+lambda*tau(1))
        call lu_solve(ID-a, Y)

        ! 5. compute K and L
        K = ky*Y
        L = ly*Y

        ! 6. Compute company (sectoral) tax revenues
        T(1) = tau(1)*q(1)*Y(1) - sum(tau*q*a(:,1))*Y(1)
        T(2) = tau(2)*q(2)*Y(2) - (1d0-lambda)*sum(tau*q*a(:,2))*Y(2)

        ! 7. check markets and budget
        markets(1) = K(1)+K(2)-Kbar
        markets(2) = p(1)*G-sum(T)-tauw*w*(Tbar-ell)-taur*r*Kbar

    end function

end module
