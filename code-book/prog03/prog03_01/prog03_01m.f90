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
    real*8, parameter :: Kbar = 10d0
    real*8, parameter :: Lbar = 20d0
    real*8, parameter :: alpha = 0.3d0
    real*8, parameter :: beta(2) = (/ 0.3d0, 0.6d0/)

contains


    ! function to determine command optimum
    function utility(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 :: utility

        ! utility function
        utility = -((x(2)**beta(1)*x(1)**(1d0-beta(1)))**alpha* &
            ((Lbar-x(2))**beta(2)*(Kbar-x(1))**(1d0-beta(2)))**(1-alpha))

    end function

end module
