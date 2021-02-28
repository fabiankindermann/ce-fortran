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

    ! number of points for interpolation
    integer, parameter :: NP = 3

    ! number of points for evaluation
    integer, parameter :: Nplot = 100

    ! lower bound for prices
    real*8, parameter :: p_min = 0.5d0

    ! upper bound for prices
    real*8, parameter :: p_max = 12.5d0

    ! marginal costs
    real*8, parameter :: c = 0.1d0

    ! declaration of variables
    real*8 :: pa(0:NP), pr(0:NP), G(0:NP, 0:NP)
    real*8 :: pa_plot(0:Nplot), pr_plot(0:Nplot)
    real*8 :: coeff_G(NP+3, NP+3), G_plot(0:Nplot, 0:Nplot)

contains

    ! the first order condition
    function foc(x_in)

        implicit none

        real*8, intent(in) :: x_in(:)
        real*8 :: foc(size(x_in, 1))

        foc(1) = 10d0 - 2d0*x_in(1) - 0.5d0*x_in(2) + 1.5d0*c
        foc(2) = 20d0 - 2d0*x_in(2) - 0.5d0*x_in(1) + c

    end function

    ! the profit function
    function profit(x_in)

        implicit none

        real*8, intent(in) :: x_in(:)
        real*8 :: profit

        profit = -(x_in(1)*(10d0 - x_in(1)) + x_in(2)*(20d0 - x_in(2) - 0.5d0*x_in(1)) &
                   - c*(10d0 - x_in(1) + 20d0 - x_in(2) - 0.5d0*x_in(1)))

    end function

end module
