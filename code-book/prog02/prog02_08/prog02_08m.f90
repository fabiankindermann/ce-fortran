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

    ! variable declaration
    implicit none
    real*8 :: eta = 1.6d0
    real*8 :: c(2) = (/0.6d0, 0.8d0/)

contains


    ! function that defines the oligopoly equations
    function cournot(q)

        implicit none
        real*8, intent(in) :: q(:)
        real*8 :: cournot(size(q, 1)),QQ
        integer :: i

        QQ = sum(q)
        do i = 1,size(q, 1)
            cournot(i) = QQ**(-1d0/eta)-1d0/eta*QQ**(-1d0/eta-1)*q(i)-c(i)*q(i)
        enddo

    end function cournot

end module
