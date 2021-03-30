!##############################################################################
! MODULE globals
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
module globals

    use toolbox
    implicit none
    integer, parameter           :: TT = 16          ! maximum age of individuals
    integer, parameter           :: NN = 1000        ! number of insured people
    integer, parameter           :: x = 8, d = 9     ! initial age and contract duration
    real*8, parameter            :: r = 0.16d0       ! interest rate
    real*8, parameter            :: xi_2 = 0.05d0    ! discount factor for premium calculation

    real*8 :: k(0:(TT-x+1)), rand, q(TT)
    real*8 :: agg_LI(NN), ind_LI(NN)
    integer :: in

contains

    subroutine init_borders()

        implicit none
        real*8 :: PP
        integer :: ij

        ! Compute right borders of intervals
        k(0) = 0d0
        PP = 1.0d0
        do ij = 1, TT-x+1
            k(ij) = k(ij-1) + q(x+ij-1)*PP
            PP = PP*(1.0d0-q(x+ij-1))
        enddo

    end subroutine

    function get_LI(rand)

        implicit none
        real*8, intent(in) :: rand
        real*8 :: get_LI
        integer :: ij

        get_LI = 0d0
        do ij = 1, TT-x+1
            if(rand < k(ij) .and. ij <= d)then
                get_LI = 1d0/(1d0+r)**ij
                exit
            endif
        enddo

    end function get_LI

end module globals
