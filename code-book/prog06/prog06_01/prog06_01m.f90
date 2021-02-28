!##############################################################################
! MODULE globals
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

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: beta = 0.9d0
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: delta = 0.0d0
    real*8, parameter :: by = 0.0d0
    real*8, parameter :: kappa = 0.0d0
    real*8, parameter :: n_p = 0.2d0
    integer, parameter :: tax = 1
    real*8 :: g(3) = (/ 0.12d0, 0.12d0, 0.0d0 /)

    ! model variables
    real*8 :: w, r, wn, Rn, p, tauw, taur, tauc, taup, pen
    real*8 :: KK, LL, YY, AA, CC, BB, GG, II
    real*8 :: a(3), c(3), util

contains


    ! equations to determine a long-run equilibrium
    function eqns(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: eqns(size(x, 1))
        real*8 :: PVI, PSI

        ! copy values to respective variables
        KK = x(1)
        if(tax == 1)tauc = x(2)
        if(tax == 2)then
            tauw = x(2)
            taur = x(2)
        endif
        if(tax == 3)tauw = x(2)
        if(tax > 3)taur = x(2)

        ! factor prices and pension payments
        r = alpha*(KK/LL)**(alpha-1d0)-delta
        w = (1d0-alpha)*(KK/LL)**alpha
        wn = w*(1d0-tauw-taup)
        Rn = 1d0+r*(1d0-taur)
        p = 1d0+tauc
        pen = kappa*w

        ! individual decisions
        PVI = wn + wn/Rn + pen/Rn**2
        PSI = p*(1d0 + (beta*Rn)**gamma/Rn &
            + (beta*Rn)**(2*gamma)/Rn**2)
        c(1) = PVI/PSI
        c(2) = (beta*Rn)**gamma*c(1)
        c(3) = (beta*Rn)**gamma*c(2)
        a(1) = 0d0
        a(2) = wn - p*c(1)
        a(3) = wn + Rn*a(2) - p*c(2)

        ! quantities
        YY = KK**alpha * LL**(1d0-alpha)
        CC = c(1) + c(2)/(1d0+n_p) + c(3)/(1d0+n_p)**2
        GG = g(1) + g(2)/(1d0+n_p) + g(3)/(1d0+n_p)**2
        AA = a(2)/(1d0+n_p) + a(3)/(1d0+n_p)**2
        II = (n_p+delta)*KK
        BB = by*YY

        ! get equations defining general equilibrium
        eqns(1) = KK + BB - AA
        eqns(2) = tauc*CC + tauw*w*LL + taur*r*AA - (r-n_p)*BB - GG

    end function

end module
