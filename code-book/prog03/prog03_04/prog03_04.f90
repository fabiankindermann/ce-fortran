!##############################################################################
! PROGRAM market3
!
! ## The static GE model with government activity
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
include "prog03_04m.f90"

program market3

    use globals
    use toolbox

    implicit none
    real*8 :: x(4), L(2), K(2), Y(2), Xd(2), Ybarn
    real*8 :: q(2), p(2), w, wn, r, rn, U, ell
    logical :: check

    ! initial guess
    x(:) = (/0.5d0, 0.5d0, 0.5d0, 0.5d0/)

    ! find market equilibrium
    call fzero(x, markets, check)

    ! check whether fzero converged
    if(check)then
        write(*,'(a/)')'Error in fzero !!!'
        stop
    endif

    ! copy producer prices
    q(1)    = 1d0
    q(2)    = x(1)
    w       = x(2)
    r       = x(3)

    ! set tax rates (uncomment respective line for different tables)
    taur    = x(4)
    !tauw    = x(4)
    !tauc(1) = x(4)
    !tauc(2) = x(4)

    ! calculate consumer prices
    p    = q*(1d0+tauc)
    wn   = w*(1d0-tauw)
    rn   = r*(1d0-taur)

    ! calculate other economic variables
    Ybarn = wn*Tbar+rn*Kbar
    Xd    = alpha*Ybarn/p
    Y(1)  = Xd(1)+G
    Y(2)  = Xd(2)
    ell   = (1d0-alpha(1)-alpha(2))*Ybarn/wn
    L     = beta*q*Y/w
    K     = (1d0-beta)*q*Y/r
    U     = Xd(1)**alpha(1)*Xd(2)**alpha(2)*ell**(1d0-alpha(1)-alpha(2))

    ! output
    write(*,'(/a)')'GOODS MARKET 1 :'
    write(*,'(3(a,f6.2,2x))')' X1 =',Xd(1),' G  =',G,' Y1 =',Y(1)
    write(*,'(3(a,f6.2,2x))')' q1 =',q(1),' p1 =',p(1),' tc1=',tauc(1)

    write(*,'(/a)')'GOODS MARKET 2 :'
    write(*,'(3(a,f6.2,2x))')' X2 =',Xd(2),' G  =',0d0,' Y2 =',Y(2)
    write(*,'(3(a,f6.2,2x))')' q2 =',q(2),' p2 =',p(2),' tc2=',tauc(2)

    write(*,'(/a)')'LABOR MARKET :'
    write(*,'(3(a,f6.2,2x))')' L1 =',L(1),' L2 =',L(2),' T-F=',Tbar-ell
    write(*,'(3(a,f6.2,2x))')' w  =',w,' wn =',wn, ' tw =',tauw

    write(*,'(/a)')'CAPITAL MARKET :'
    write(*,'(3(a,f6.2,2x))')' K1 =',K(1),' K2 =',K(2),' K  =',Kbar
    write(*,'(3(a,f6.2,2x))')' r  =',r,' rn =',rn, ' tr =',taur

    write(*,'(/a)')'GOVERNMENT :'
    write(*,'(6(a,f6.2,2x))')' tc1=',tauc(1)*q(1)*Xd(1), &
        ' tc2=',tauc(2)*q(2)*Xd(2),' tw =',tauw*w*(Tbar-ell), &
        ' tr =',taur*r*Kbar,' G  =',q(1)*G

    write(*,'(/a)')'UTILITY :'
    write(*,'(a,f6.2,2x)')' U  =',U

end program
