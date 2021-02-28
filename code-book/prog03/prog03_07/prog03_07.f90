!##############################################################################
! PROGRAM market6
!
! ## The static GE model with international trade
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
include "prog03_07m.f90"

program market6

    use globals
    use toolbox

    implicit none
    real*8 :: x(4), U(2)
    integer :: i
    logical :: check

    ! initial guess
    x(:) = (/1d0, 0.5d0, 0d0, 0d0/)

    ! find market equilibrium
    call fzero(x, markets, check)

    ! check whether fzero converged
    if(check)then
        write(*,'(a/)')'Error in fzero !!!'
        stop
    endif

    do i = 1, 2

        ! get utility level
        U(i) = Xd(1,i)**alpha(1)*Xd(2,i)**alpha(2)*ell(i)**(1d0-alpha(1)-alpha(2))

        if(i == 1)then
            write(*,'(/a)')'COUNTRY A:'
        else
            write(*,'(a)')'COUNTRY B:'
        endif

        ! output
        write(*,'(/a)')'GOODS MARKET 1 :'
        write(*,'(5(a,f6.2,2x))')' X11=',a(1,1)*Y(1,i),' X12=',a(1,2)*Y(2,i), &
            ' X1 =',Xd(1,i),' G  =',G,' Y1 =',Y(1,i)
        write(*,'(4(a,f6.2,2x))')' q1 =',q(1,i),' p1 =',p(1,i),' tc1=',tauc(1,i)

        write(*,'(/a)')'GOODS MARKET 2 :'
        write(*,'(6(a,f6.2,2x))')' X21=',a(2,1)*Y(1,i),' X22=',a(2,2)*Y(2,i), &
            ' X2 =',Xd(2,i),' G  =',0d0,' Y2 =',Y(2,i),&
            ' EXP=',Y(2,i)-Xd(2,i)-a(2,1)*Y(1,i)-a(2,2)*Y(2,i)
        write(*,'(4(a,f6.2,2x))')' q2 =',q(2,i),' p2 =',p(2,i),' tc2=',tauc(2,i)

        write(*,'(/a)')'LABOR MARKET :'
        write(*,'(3(a,f6.2,2x))')' L1 =',L(1,i),' L2 =',L(2,i),' T-l=',Tbar(i)-ell(i)
        write(*,'(3(a,f6.2,2x))')' w  =',w(i),' wn =',wn(i), ' tw =',tauw(i)

        write(*,'(/a)')'CAPITAL MARKET :'
        write(*,'(4(a,f6.2,2x))')' K1 =',K(1,i),' K2 =',K(2,i),' K  =',Kbar(i), &
            ' EXP=',Kbar(i)-K(1,i)-K(2,i)
        write(*,'(3(a,f6.2,2x))')' r  =',r,' rn =',rn(i), ' tr =',taur(i)

        write(*,'(/a)')'GOVERNMENT :'
        write(*,'(6(a,f6.2,2x))')' tc1=',tauc(1,i)*q(1,i)*Xd(1,i), &
            ' tc2=',tauc(2,i)*q(2,i)*Xd(2,i),' tw =',tauw(i)*w(i)*(Tbar(i)-ell(i)), &
            ' tr =',taur(i)*r*Kbar(i),' G  =',q(1,i)*G

        write(*,'(/a)')'UTILITY :'
        write(*,'(a,f6.2,2x///)')' U  =',U(i)

        ! IO-Table
        write(*,'(a)')' ---------------------------------------------'
        write(*,'(a,f6.2,2x,f6.2,a,3(f6.2,2x),a,f6.2)')' | ',q(1,i)*a(1,1)*Y(1,i), &
            q(1,i)*a(1,2)*Y(2,i),' | ',q(1,i)*Xd(1,i),q(1,i)*G, &
            q(1,i)*(Y(1,i)-Xd(1,i)-a(1,1)*Y(1,i)-a(1,2)*Y(2,i)-G),' | ',q(1,i)*Y(1,i)
        write(*,'(a,f6.2,2x,f6.2,a,3(f6.2,2x),a,f6.2)')' | ',q(2,i)*a(2,1)*Y(1,i), &
            q(2,i)*a(2,2)*Y(2,i),' | ',q(2,i)*Xd(2,i),0d0, &
            q(2,i)*(Y(2,i)-Xd(2,i)-a(2,1)*Y(1,i)-a(2,2)*Y(2,i)),' | ',q(2,i)*Y(2,i)
        write(*,'(a)')' ---------------------------------------------'
        write(*,'(a,f6.2,2x,f6.2,a)')' | ',w(i)*L(1,i),w(i)*L(2,i),' | '
        write(*,'(a,f6.2,2x,f6.2,a)')' | ',r*K(1,i),r*K(2,i),' | '
        write(*,'(a)')' ------------------'
        write(*,'(a,f6.2,2x,f6.2//)')'   ',q(1,i)*Y(1,i),q(2,i)*Y(2,i)
    enddo

end program
