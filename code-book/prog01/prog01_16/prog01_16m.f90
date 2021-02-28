!##############################################################################
! MODULE Globals
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
module Globals

    implicit none

    ! time preference
    real*8 :: beta

    ! risk aversion
    real*8 :: eta

end module


!##############################################################################
! MODULE UtilFunc
!##############################################################################
module UtilFunc

    implicit none

contains


    ! a utility function
    function utility(c1, c2)

        use Globals

        implicit none
        real*8, intent(in) :: c1, c2
        real*8 :: utility

        utility = 1d0/(1d0-eta)*c1**(1d0-eta) + beta/(1d0-eta)*c2**(1d0-eta)

    end function

end module
