!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    implicit none
    real*8, parameter :: Kbar = 10d0
    real*8, parameter :: Tbar = 30d0
    real*8, parameter :: alpha(2) = (/0.3d0, 0.4d0/)
    real*8, parameter :: beta(2) = (/ 0.3d0, 0.6d0/)
    real*8, parameter :: gamma(2) = (/ -2.5d0, 0.5d0/)
    real*8, parameter :: G = 3d0
    real*8, parameter :: B = 0.1d0
    real*8            :: tauw = 0d0
    real*8            :: tauu = 0d0
    real*8            :: taur = 0d0
    real*8            :: tauc(2) = 0d0

contains


    ! function to determine market equilibrium
    function markets(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: markets(size(x, 1))
        real*8 :: Ybarn, q(2), p(2), PX, w, wn, r, rn, ur

        ! copy producer prices and taxes
        q(1)    = 1d0
        q(2)    = x(1)
        w       = x(2)
        r       = x(3)

        ! set tax rates (uncomment respective line for different tables)
        taur    = x(4)
        !tauw    = x(4)
        !tauc(1) = x(4)
        !tauc(2) = x(4)

        ! copy the unempoyment rate
        ur      = x(5)

        ! calculate consumer prices and total income
        p     = q*(1d0+tauc)
        PX    = p(1)*alpha(1)/sum(alpha)+p(2)*alpha(2)/sum(alpha)
        tauu  = PX*B*ur/(w*(1-ur))
        wn    = w*(1d0-tauw-tauu)*(1-ur)+PX*B*ur
        rn    = r*(1d0-taur)
        Ybarn = wn*Tbar+rn*Kbar

        ! get market equations
        markets(1) = alpha(1)*Ybarn/p(1)+G-(beta(1)/w)**beta(1)* &
            ((1d0-beta(1))/r)**(1d0-beta(1))*q(1)*(alpha(1)*Ybarn/p(1)+G)
        markets(2) = 1d0-(beta(2)/w)**beta(2)*((1d0-beta(2))/r)**(1d0-beta(2))*q(2)
        markets(3) = beta(1)/w*q(1)*(alpha(1)*Ybarn/p(1)+G)+ &
            beta(2)/w*q(2)*alpha(2)*Ybarn/p(2)+(1-ur)*((1d0-alpha(1)-alpha(2))*Ybarn/wn-Tbar)
        markets(4) = q(1)*G-tauc(1)/(1d0+tauc(1))*alpha(1)*Ybarn- &
            tauc(2)/(1d0+tauc(2))*alpha(2)*Ybarn-&
            tauw*w*(1-ur)*(Tbar-(1d0-alpha(1)-alpha(2))/wn*Ybarn)-taur*r*Kbar
        markets(5) = log(w*(1-tauw-tauu)/PX)-gamma(1)+gamma(2)*log(ur)

    end function markets

end module

