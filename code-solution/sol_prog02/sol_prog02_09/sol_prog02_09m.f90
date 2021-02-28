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

    ! model parameters
    real*8, parameter :: alpha = 1d0
    real*8, parameter :: eta = 1.5d0
    integer, parameter :: m = 3

    ! numerical parameters
    real*8, parameter :: p_l = 0.1d0
    real*8, parameter :: p_u = 3.0d0

    ! variables
    integer, parameter :: N = 10
    integer, parameter :: NP = 1000
    real*8 :: p_equi, P(0:N), q(0:N), D(0:N)
    real*8 :: coeff_q(N + 3), coeff_d(N + 3), p_plot(0:NP), q_s_plot(0:NP), q_d_plot(0:NP)

    ! communicatation variables
    integer :: ip_com

contains

    ! the first order condition
    function foc(q_in)

        real*8, intent(in) :: q_in
        real*8 :: foc

        foc = p(ip_com) - q_in*p(ip_com)**(1d0+eta)/eta - alpha*sqrt(q_in) - q_in**2d0

    end function

    ! the market equilibrium
    function market(p_in)

        use toolbox

        real*8, intent(in) :: p_in
        real*8 :: market

        market = p_in**(-eta) - m*spline_eval(p_in, coeff_q, p_l, p_u)

    end function

end module
