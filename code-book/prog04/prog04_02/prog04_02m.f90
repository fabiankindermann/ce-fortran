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
    integer, parameter :: TT = 5        ! number of time periods
    integer, parameter :: NN = 3        ! number of stocks
    real*8, parameter :: gamma = 10d0   ! risk aversion
    real*8, parameter :: r_f = 0.05     ! the risk free rate

    real*8 :: mu(NN), sig(NN,NN)

contains


    ! the first order conditions of the investor
    function focs(x_in)

        implicit none
        real*8, intent(in) :: x_in(:)
        real*8 :: focs(size(x_in, 1))
        real*8 :: omega(NN), omega_f, eta(NN), eta_f

        ! copy portfolio weights
        omega = x_in(1:NN)
        omega_f = 1d0 - sum(omega)

        ! copy lagrangean multipliers
        eta   = x_in(NN+1:2*NN)
        eta_f = x_in(2*NN+1)

        ! set up first order conditions
        focs(1:NN)      = mu - r_f - gamma*matmul(sig, omega) + eta - eta_f
        focs(NN+1:2*NN) = omega + eta - sqrt(omega**2 + eta**2)
        focs(2*NN+1)    = omega_f + eta_f - sqrt(omega_f**2 + eta_f**2)

    end function

end module

