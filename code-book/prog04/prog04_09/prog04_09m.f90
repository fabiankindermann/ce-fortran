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
    integer, parameter :: XX = 16            ! maximum age of individuals
    integer, parameter :: TT = XX - 1        ! maximum number of periods
    integer, parameter :: KK = 100000        ! number of mortality path simulations
    integer, parameter :: NN = 5             ! number of insurance contracts

    real*8, parameter  :: sigma_eps = 2.0d0  ! variance of the innovation to the random walk
    real*8, parameter  :: c = -1.0d0         ! drift of random walk

    real*8, parameter  :: r = 0.16d0         ! interest rate
    real*8, parameter  :: xi_1 = 0.05d0      ! loading factor for premium calculation (annuity)
    real*8, parameter  :: xi_2 = 0.05d0      ! discount factor for premium calculation (life insurance)
    real*8, parameter  :: gamma = 10.0d0     ! risk aversion

    real*8 :: q(XX, 0:TT, 0:KK), alpha(XX), beta(XX), kappa(0:TT), eps(TT)
    real*8 :: PP(XX, 0:TT, 0:KK), LE(XX, 0:KK)
    real*8 :: p_m(NN), omega(NN), n_new(NN)
    real*8 :: l(NN, KK), l_p(KK)
    real*8 :: mu(NN), sig(NN, NN), mu_p, sig_p, sum_premia
    real*8 :: VaR, CVaR, IID(NN), x_in(2*NN+1)
    real*8 :: age(0:TT+1), PP_av(0:TT+1), PP_min(0:TT+1), PP_max(0:TT+1)
    integer :: ix, ik, it, ik_max, ik_min, ii, im, id, ij, ii_VaR
    logical :: check

    ! define an insurance type
    type :: insurance_type
        character(len=2) :: itype
        integer :: ix, im, id
        real*8 :: BB, n_ins
    end type

    ! declare the insurance variable
    type(insurance_type) :: ins(NN)

contains


    ! first order condition of the insurer
    function focs(x_in)

        implicit none
        real*8, intent(in) :: x_in(:)
        real*8 :: focs(size(x_in, 1))
        real*8 :: omega(NN), eta(NN), lamda

        ! copy portfolio weights
        omega = x_in(1:NN)

        ! copy lagrangean multipliers
        lamda = x_in(NN+1)
        eta   = x_in(NN+2:2*NN+1)

        ! set up first order conditions
        focs(1:NN) = - mu - gamma*matmul(sig, omega) - lamda + eta
        focs(NN+1) = sum(omega) - 1d0
        focs(NN+2:2*NN+1) = omega + eta - sqrt(omega**2 + eta**2)

    end function focs

end module globals


