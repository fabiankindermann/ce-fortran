!##############################################################################
! MODULE globals
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
module globals

    implicit none

    ! declaration of variables
    integer, parameter :: n = 11
    real*8, parameter :: gamma = 0.5d0

contains

    ! subroutine to compute utilities within a certain interval
    subroutine utility_int(a_local, b_local, u_local)

        implicit none

        ! input variables
        real*8, intent(in) :: a_local, b_local
        real*8, intent(out) :: u_local(:)
        real*8 :: c_local(size(u_local))
        integer :: n, j

        ! check whether 0 <= a < b:
        if(a_local <= 0d0 .or. b_local < a_local)then
            stop 'Error: a <= 0 or b < a: subroutine can`t be used'
        endif

        ! size of the array u_local
        n = size(u_local)

        ! calculation of consumption level and utilities
        do j = 1, n
            c_local(j) = a_local + dble(j-1)/dble(n-1)*(b_local - a_local)
            u_local(j) = utility(c_local(j))
        enddo

    end subroutine

    ! subroutine to compute utility
    function utility(c_local)

        implicit none

        ! input variables
        real*8, intent(in) :: c_local
        real*8 :: utility

        ! check whether c_local <= 0; otherwise print an error message
        if(c_local <= 0d0)then
            stop 'Error: c_local <= 0 in function utility'
        endif

        ! compute utility
        utility = c_local**(1d0-1d0/gamma)/(1d0-1d0/gamma)

    end function

end module
