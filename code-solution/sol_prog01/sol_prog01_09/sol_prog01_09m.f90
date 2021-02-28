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

    ! number of gridpoints for consumption
    integer, parameter :: NC = 100

    ! number of different values for gamma
    integer, parameter :: NG = 4

    ! initialize array for different vaues of gamma
    real*8, parameter :: gamma_array(NG) = (/0.25d0, 0.5d0, 0.75d0, 1.25d0/)

contains

    ! function to compute utility
    function utility(c_local, gamma_local)

        ! input variables
        real*8, intent(in) :: c_local, gamma_local
        real*8 :: utility

        ! check whether c_loc > 0; otherwise print an error message
        if(c_local <= 0d0)then
            stop 'Error: c <= 0 in function utility'
        endif

        ! compute utility
        utility = c_local**(1d0 - 1d0/gamma_local)/(1d0 - 1d0/gamma_local)

    end function

end module
