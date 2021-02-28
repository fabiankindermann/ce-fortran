!##############################################################################
! MODULE GLOBALS
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
module globals

    implicit none
    real*8, parameter :: Kbar = 10d0
    real*8, parameter :: Tbar = 30d0
    real*8, parameter :: alpha(2) = (/0.3d0, 0.4d0/)
    real*8, parameter :: beta(2) = (/ 0.3d0, 0.6d0/)
    real*8, parameter :: a0(2) = (/0.2d0, 0.2d0/)
    real*8, parameter :: a(2, 2) = reshape((/ 0d0, 0.3d0, 0.2d0, 0d0 /),(/2,2/))
    real*8, parameter :: ID(2, 2) = reshape((/ 1d0, 0d0, 0d0, 1d0 /),(/2,2/))
    real*8, parameter :: G = 3d0
    real*8            :: tauw = 0d0
    real*8            :: taur = 0d0
    real*8            :: tauc(2) = 0d0
    real*8            :: Ybarn, q(2), p(2), w, wn, r, rn
    real*8            :: Xd(2), Y(2), ky(2), ly(2), K(2), L(2), ell

contains


    ! function to determine market equilibrium
    function markets(x)

        use toolbox

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: markets(size(x, 1))

        ! copy producer prices and taxes
        w       = 1d0
        r       = x(1)
        tauc(1) = x(2)
        tauc(2) = tauc(1)

        ! 1. calkulate K/Y and L/Y
        ky = a0*((1d0-beta)/beta*w/r)**beta
        ly = a0*(beta/(1d0-beta)*r/w)**(1d0-beta)

        ! 2. determine producer prices
        q = w*ly+r*ky
        call lu_solve(ID-transpose(a), q)

        ! 3. consumer prices and demands
        p  = q*(1d0+tauc)
        wn = w*(1d0-tauw)
        rn = r*(1d0-taur)
        Ybarn = wn*Tbar+rn*Kbar
        Xd = alpha/p*Ybarn
        ell  = (1d0-alpha(1)-alpha(2))/wn*Ybarn

        ! 4. determine output levels
        Y(1) = Xd(1)+G
        Y(2) = Xd(2)
        call lu_solve(ID-a, Y)

        ! 5. compute K and L
        K = ky*Y
        L = ly*Y

        ! 6. check markets and budget
        markets(1) = K(1)+K(2)-Kbar
        markets(2) = q(1)*G-sum(tauc*q*Xd)-tauw*w*(Tbar-ell)-taur*r*Kbar

    end function

end module
