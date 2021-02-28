!##############################################################################
! PROGRAM market5
!
! ## The static GE model with intermediate goods
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
include "prog03_06m.f90"

program market5

    use globals
    use toolbox

    implicit none
    real*8 :: x(2), U
    logical :: check

    ! initial guess
    x(:) = (/0.2d0, 0d0/)

    ! find market equilibrium
    call fzero(x, markets, check)

    ! check whether fzero converged
    if(check)then
        write(*,'(a/)')'Error in fzero !!!'
        stop
    endif

    ! get utility level
    U = Xd(1)**alpha(1)*Xd(2)**alpha(2)*ell**(1d0-alpha(1)-alpha(2))

    ! output
    write(*,'(/a)')'GOODS MARKET 1 :'
    write(*,'(5(a,f6.2,2x))')' X11=',a(1,1)*Y(1),' X12=',a(1,2)*Y(2), &
        ' X1 =',Xd(1),' G  =',G,' Y1 =',Y(1)
    write(*,'(4(a,f6.2,2x))')' q1 =',q(1),' p1 =',p(1),' tc1=',tauc(1)

    write(*,'(/a)')'GOODS MARKET 2 :'
    write(*,'(5(a,f6.2,2x))')' X21=',a(2,1)*Y(1),' X22=',a(2,2)*Y(2), &
        ' X2 =',Xd(2),' G  =',0d0,' Y2 =',Y(2)
    write(*,'(4(a,f6.2,2x))')' q2 =',q(2),' p2 =',p(2),' tc2=',tauc(2)

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
    write(*,'(a,f6.2,2x///)')' U  =',U

    ! IO-Table
    write(*,'(a)')'IO-TABLE:'
    write(*,'(a)')' -----------------------------------'
    write(*,'(a,2(f6.2,2x,f6.2,a),f6.2)')' | ',q(1)*a(1,1)*Y(1), &
        q(1)*a(1,2)*Y(2),' | ',q(1)*Xd(1),q(1)*G,' | ',q(1)*Y(1)
    write(*,'(a,2(f6.2,2x,f6.2,a),f6.2)')' | ',q(2)*a(2,1)*Y(1), &
        q(2)*a(2,2)*Y(2),' | ',q(2)*Xd(2),0d0,' | ',q(2)*Y(2)
    write(*,'(a)')' -----------------------------------'
    write(*,'(a,f6.2,2x,f6.2,a)')' | ',w*L(1),w*L(2),' | '
    write(*,'(a,f6.2,2x,f6.2,a)')' | ',r*K(1),r*K(2),' | '
    write(*,'(a)')' ------------------'
    write(*,'(a,f6.2,2x,f6.2)')'   ',q(1)*Y(1),q(2)*Y(2)

end program
