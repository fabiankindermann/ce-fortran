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

    ! commons
    real*8, parameter :: alpha(2) = (/0.3d0, 0.4d0/)
    real*8, parameter :: a0(2) = (/0.2d0, 0.2d0/)
    real*8, parameter :: a(2, 2) = reshape((/ 0d0, 0.3d0, 0.2d0, 0d0 /),(/2,2/))
    real*8, parameter :: G = 3d0
    real*8, parameter :: ID(2, 2) = reshape((/ 1d0, 0d0, 0d0, 1d0 /),(/2,2/))
    real*8            :: r

    ! country specific
    real*8, parameter :: Tbar(2) = (/30d0, 10d0/)
    real*8, parameter :: Kbar(2) = (/10d0, 30d0/)
    real*8, parameter :: beta(2, 2) = reshape((/ 0.3d0, 0.6d0, 0.3d0, 0.6d0 /),(/2,2/))
    real*8            :: tauw(2) = (/0d0, 0d0/)
    real*8            :: taur(2) = (/0d0, 0d0/)
    real*8            :: tauc(2,2) = 0d0
    real*8            :: Ybarn(2), q(2,2), p(2,2), w(2), wn(2), rn(2)
    real*8            :: Xd(2,2), Y(2,2), ky(2,2), ly(2,2), K(2,2), L(2,2), ell(2)

contains


    ! function to determine market equilibrium
    function markets(x)

        use toolbox

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: markets(size(x, 1))
        real*8 :: mat(4, 4), vec(4)
        integer :: i

        ! copy producer prices and taxes
        w(1)      = 1d0
        w(2)      = x(1)
        r         = x(2)
        tauc(1,1) = x(3)
        tauc(1,2) = x(4)
        tauc(2,:) = tauc(1,:)

        ! 1. calculate K/Y and L/Y
        do i = 1, 2
            ky(:,i) = a0*((1d0-beta(:,i))/beta(:,i)*w(i)/r)**beta(:,i)
            ly(:,i) = a0*(beta(:,i)/(1d0-beta(:,i))*r/w(i))**(1d0-beta(:,i))
        enddo

        ! 2. determine producer prices
        q(:,1) = w(1)*ly(:,1)+r*ky(:,1)
        call lu_solve(ID-transpose(a), q(:,1))
        q(1,2) = (a(2,1)*q(2,1)+w(2)*ly(1,2)+r*ky(1,2))/(1d0-a(1,1))
        q(2,2) = q(2,1)

        ! 3. consumer prices and demands
        do i = 1, 2
            p(:,i) = q(:,i)*(1d0+tauc(:,i))
            wn(i) = w(i)*(1d0-tauw(i))
            rn(i) = r*(1d0-taur(i))
            Ybarn(i) = wn(i)*Tbar(i)+rn(i)*Kbar(i)
            Xd(:,i) = alpha/p(:,i)*Ybarn(i)
            ell(i)  = (1d0-alpha(1)-alpha(2))/wn(i)*Ybarn(i)
        enddo

        ! 4. determine output levels
        vec(1) = Xd(1,1)+G
        vec(2) = Xd(2,1)+Xd(2,2)
        vec(3) = Xd(1,2)+G
        vec(4) = Tbar(1)-ell(1)

        mat(1, :) = (/1d0-a(1,1), -a(1,2), 0d0, 0d0/)
        mat(2, :) = (/-a(2,1), 1d0-a(2,2), -a(2,1), 1d0-a(2,2)/)
        mat(3, :) = (/0d0, 0d0, 1d0-a(1,1), -a(1,2)/)
        mat(4, :) = (/ly(1,1), ly(2,1), 0d0, 0d0/)

        call lu_solve(mat, vec)
        Y(1,1) = vec(1)
        Y(2,1) = vec(2)
        Y(1,2) = vec(3)
        Y(2,2) = vec(4)

        ! 5. compute K and L
        K = ky*Y
        L = ly*Y

        ! 6. check markets and budget
        markets(1) = L(1,2)+L(2,2)-(Tbar(2)-ell(2))
        markets(2) = sum(K)-sum(Kbar)
        markets(3) = q(1,1)*G-sum(tauc(:,1)*q(:,1)*Xd(:,1))-&
            tauw(1)*w(1)*(Tbar(1)-ell(1))-taur(1)*r*Kbar(1)
        markets(4) = q(1,2)*G-sum(tauc(:,2)*q(:,2)*Xd(:,2))-&
            tauw(2)*w(2)*(Tbar(2)-ell(2))-taur(2)*r*Kbar(2)

    end function

end module

