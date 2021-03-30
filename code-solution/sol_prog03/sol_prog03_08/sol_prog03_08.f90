!##############################################################################
! PROGRAM factortax
!
! ## The static GE model with sector-specific factor input tax
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog03_08m.f90"

program factortax

    use globals
    use toolbox

    implicit none
    real*8 :: x(2)
    logical :: check

    ! initial guess
    x(:) = (/0.2d0, 0d0/)
    !tauk(1) = 0.2d0
    taul(2) = 0.2d0

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
    write(*,'(5(a,f6.2,2x))')' q1 =',q(1),' p1 =',p(1),' tc1=',tauc(1), &
        ' tl =',taul(1), ' tk =',tauk(1)

    write(*,'(/a)')'GOODS MARKET 2 :'
    write(*,'(5(a,f6.2,2x))')' X21=',a(2,1)*Y(1),' X22=',a(2,2)*Y(2), &
        ' X2 =',Xd(2),' G  =',0d0,' Y2 =',Y(2)
    write(*,'(5(a,f6.2,2x))')' q2 =',q(2),' p2 =',p(2),' tc2=',tauc(2), &
        ' tl =',taul(2), ' tk =',tauk(2)

    write(*,'(/a)')'LABOR MARKET :'
    write(*,'(3(a,f6.2,2x))')' L1 =',L(1),' L2 =',L(2),' T-l=',Tbar-ell
    write(*,'(4(a,f6.2,2x))')' w  =',w,' wn =',wn, ' tw =',tauw, ' DIFF =', &
            L(1)+L(2)+ell-Tbar

    write(*,'(/a)')'CAPITAL MARKET :'
    write(*,'(3(a,f6.2,2x))')' K1 =',K(1),' K2 =',K(2),' K  =',Kbar
    write(*,'(4(a,f6.2,2x))')' r  =',r,' rn =',rn, ' tr =',taur, ' DIFF =', &
            K(1)+K(2)-Kbar

    write(*,'(/a)')'GOVERNMENT :'
    write(*,'(7(a,f6.2,2x))')' tc1=',tauc(1)*q(1)*Xd(1), &
        ' tc2=',tauc(2)*q(2)*Xd(2),' tw =',tauw*w*(Tbar-ell), &
        ' tr =',taur*r*Kbar,' tl =',sum(taul*w*L),' tk =',sum(tauk*r*K)
    write(*,'(2(a,f6.2,2x))') ' G  =',q(1)*G, ' DIFF =', q(1)*G-sum(tauc*q*Xd)- &
            tauw*w*(Tbar-ell)-taur*r*Kbar-sum(taul*w*L)-sum(tauk*r*K)

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
    write(*,'(a,f6.2,2x,f6.2,a)')' | ',w*(1+taul(1))*L(1),w*(1+taul(2))*L(2),' | '
    write(*,'(a,f6.2,2x,f6.2,a)')' | ',r*(1+tauk(1))*K(1),r*(1+tauk(2))*K(2),' | '
    write(*,'(a)')' ------------------'
    write(*,'(a,f6.2,2x,f6.2)')'   ',q(1)*Y(1),q(2)*Y(2)


end program
