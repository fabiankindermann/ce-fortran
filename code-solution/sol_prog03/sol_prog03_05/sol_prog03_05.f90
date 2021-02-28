!##############################################################################
! PROGRAM marketCESTAX
!
! ## The static general equilibrium model with CES production and Government
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
include "sol_prog03_05m.f90"

program marketCESTAX

    use globals
    use toolbox

    implicit none
    real*8 :: x(2)
    logical :: check

    ! initial guess
    x(:) = (/0.5d0, 0.4d0/)

    ! find market equilibrium
    call fzero(x, markets, check)

    ! check whether fzero converged
    if(check)then
        write(*,'(a/)')'Error in fzero !!!'
        stop
    endif

    ! get utility level
    U = ((1d0-alpha)**(1d0/nu)*C**mu+alpha**(1d0/nu)*ell**mu)**(1d0/mu)

    ! output
    write(*,'(/a)')'GOODS MARKET 1 :'
    write(*,'(4(a,f6.2,2x))')' X1 =',Y(1),' Y1 =',Y(1), &
        ' q1 =',p(1),' p1 =',p(1)

    write(*,'(/a)')'GOODS MARKET 2 :'
    write(*,'(4(a,f6.2,2x))')' X2 =',Y(2),' Y2 =',Y(2), &
        ' q2 =',p(2),' p2 =',p(2)

    write(*,'(/a)')'LABOR MARKET :'
    write(*,'(4(a,f6.2,2x))')' L1 =',L(1),' L2 =',L(2),' L  =',Tbar-ell, &
        ' w  =',w

    write(*,'(/a)')'CAPITAL MARKET :'
    write(*,'(4(a,f6.2,2x))')' K1 =',K(1),' K2 =',K(2),' K  =',Kbar, &
        ' r  =',r

    write(*,'(/a)')'GOVERNMENT :'
    write(*,'(6(a,f6.2,2x))')' tc1=',tauc(1)*q(1)*Xd(1), &
        ' tc2=',tauc(2)*q(2)*Xd(2),' tw =',tauw*w*(Tbar-ell), &
        ' tr =',taur*r*Kbar,' G  =',q(1)*G

    write(*,'(/a)')'UTILITY :'
    write(*,'(a,f6.2,2x)')' U  =',U

end program
