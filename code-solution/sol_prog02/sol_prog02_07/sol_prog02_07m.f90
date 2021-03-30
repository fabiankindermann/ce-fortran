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

    implicit none

    ! number of nodes for interpolation
    integer, parameter :: n = 2

    ! grid to evaluate plot
    integer, parameter :: n_plot = 10

    ! lower bound for interpolation
    real*8, parameter :: tau_l = 35d0

    ! upper bound for interpolation
    real*8, parameter :: tau_u = 45d0

    ! declaration of variables
    real*8 ::  tau_data(0:n), t_data(0:n)
    real*8 :: tau_plot(0:n_plot), t_plot(0:n_plot)

contains

    ! the tax function that should be maximized
    function tax_func(tau_in)

        use toolbox

        implicit none

        ! declaration of variables
        real*8, intent(in) :: tau_in
        real*8 :: tax_func

        ! set up the tax function we want to maximize
        tax_func = -poly_interpol(tau_in, tau_data, t_data)

    end function

end module
