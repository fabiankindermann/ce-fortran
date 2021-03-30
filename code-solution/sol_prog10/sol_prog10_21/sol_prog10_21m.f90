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

    ! number of years the household retires
    integer, parameter :: JR = 45

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 5

    ! number of transitory (epsilon) shocks
    integer, parameter :: NS = 5

    ! number of rate of return (vtheta) shocks
    integer, parameter :: NR = 5

    ! number of eps-vtheta shocks
    integer, parameter :: NSR = NS*NR

    ! number of points on the asset grid
    integer, parameter :: NA = 80

    ! number of points on the cash-on-hand-grid
    integer, parameter :: NX = 80

    ! number of points on the annuity grid
    integer, parameter :: NY = NA

    ! household preference parameters
    real*8, parameter :: gamma = 0.10d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: nu = 0d0!10d0

    ! risk processes
    real*8, parameter :: sigma_zeta   = 0.0738d0
    real*8, parameter :: sigma_eps    = 0.0106d0
    real*8, parameter :: sigma_vtheta = 0.157d0**2d0
    real*8, parameter :: rho          = 0.00d0

    ! risk free rate and risk premium
    real*8, parameter :: r_f  = 0.02d0
    real*8, parameter :: mu_r = 0.04d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 250d0
    real*8, parameter :: a_grow = 0.07d0

    ! size of the retirement savings account grid
    real*8, parameter :: ya_l = 0d0
    real*8, parameter :: ya_u = 30d0
    real*8, parameter :: ya_grow = a_grow
    real*8, parameter :: xi = 0.07d0

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.05d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! lower and upper bound for X-grid
    real*8 :: X_l, X_u

    ! price of the annuity
    real*8 :: p_a(JJ)

    ! lower and upper bound for the eta-grid
    real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w, eff(JJ), pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ+1), ya_coh(JJ), o_coh(JJ+1)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_a(JJ), cv_ya(JJ), cv_o(JJ)

    ! individual variables
    real*8 :: X(0:NX), a(0:NA), ya(0:NY), c(JJ, 0:NX, 0:NY)
    real*8 :: a_plus(JJ, 0:NX, 0:NY), omega_plus(JJ, 0:NA, 0:NY), omegay_plus(JJ, 0:NA, 0:NY)
    real*8 :: V(JJ, 0:NX, 0:NY) = 0d0, Q(JJ, 0:NA, 0:NY) = 0d0, S(JJ, 0:NA, 0:NY) = 0d0
    real*8 :: phi_X(JJ, 0:NX, 0:NY), phi_a(JJ, 0:NA, 0:NY), phi_r(JJ, 0:NA, 0:NY)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, 0:NY), S_ya(JJ, 0:NA, 0:NY), V_ya(JJ, 0:NX, 0:NY) = 0d0
    real*8 :: Q_al(1:JJ, 0:NA, 0:NY), Q_ya(1:JJ, 0:NA, 0:NY)
    integer :: ij_com, ix_com, ia_com, ir_com
    real*8 :: cons_com

contains


    ! the first order condition regarding consumption
    function foc_cons(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_cons, a_plus, varphi, tomorrow
        integer :: ial, iar

        ! calculate tomorrows assets
        a_plus  = x_in

        ! calculate consumption
        cons_com = X(ix_com) - a_plus

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com, ial, ir_com) + &
                       (1d0-varphi)*RHS(ij_com, iar, ir_com)

        ! calculate first order condition for consumption
        foc_cons = cons_com - tomorrow

    end function


    ! the first order condition regarding portfolio choice
    function foc_port(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_port
        real*8 :: omega_p, R_port, X_p,  c_p, varphi_X, dist
        real*8 :: ya_p, varphi_y
        integer :: ixl, ixr, iyl, iyr, iw, isr

        ! store portfolio share
        omega_p  = x_in
        foc_port = 0d0

        if(ij_com == JJ)then
            do isr = 1, NSR

                ! get future cash-on-hand and interpolate
                R_port = 1d0 + (r_f + omega_p*(mu_r + vtheta(isr)))

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! get FOC
                foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*nu*margu(R_port*a(ia_com))
            enddo
        elseif(ij_com >= JR-1)then
            do isr = 1, NSR

                ! get future cash-on-hand and interpolate
                R_port = 1d0 + (r_f + omega_p*(mu_r + vtheta(isr)))
                X_p = R_port*a(ia_com) + pen(ij_com+1) + ya(ir_com)
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! calculate consumption and FOC
                c_p = varphi_X      *c(ij_com+1, ixl, ir_com) + &
                      (1d0-varphi_X)*c(ij_com+1, ixr, ir_com)
                c_p = max(c_p, 1d-10)
                foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*(psi(ij_com+1)*margu(c_p)+ &
                    nu*(1d0-psi(ij_com+1))*margu(R_port*a(ia_com)))
            enddo
        else
            do iw = 1, NW
                do isr = 1, NSR

                    ! get future cash-on-hand and interpolate
                    R_port = 1d0 + (r_f + omega_p*(mu_r + vtheta(isr)))
                    X_p = R_port*a(ia_com)/eps(isr) + w*eff(ij_com+1)*zeta(iw)
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                    ! get future balance of retirement account and interpolate
                    ya_p = ya(ir_com)/eps(isr)
                    call linint_Grow(ya_p, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_epsvtheta(isr)

                    ! calculate consumption and FOC
                    if(varphi_X <= varphi_y)then
                        c_p = varphi_X           *c(ij_com+1, ixl, iyl) + &
                              (varphi_y-varphi_X)*c(ij_com+1, ixr, iyl) + &
                              (1d0-varphi_y)     *c(ij_com+1, ixr, iyr)
                    else
                        c_p = varphi_y           *c(ij_com+1, ixl, iyl) + &
                              (varphi_X-varphi_y)*c(ij_com+1, ixl, iyr) + &
                              (1d0-varphi_X)     *c(ij_com+1, ixr, iyr)
                    endif
                    c_p = max(c_p, 1d-10)
                    foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*(psi(ij_com+1)*margu(eps(isr)*c_p) + &
                        nu*(1d0-psi(ij_com+1))*margu(R_port*a(ia_com)))
                enddo
            enddo
        endif

    end function


    ! the first order condition regarding retirement savings problem
    function foc_annu(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_annu
        real*8 :: omega_r, al_p, ya_p, Qal, Qya
        real*8 :: varphi_a, varphi_y
        integer :: ial, iar, iyl, iyr

        ! store share invested in retirement account
        omega_r  = x_in

        ! determine future liquid wealth and pre-return retirement wealth
        al_p = (1d0-omega_r)*a(ia_com)
        ya_p = ya(ir_com) + omega_r*a(ia_com)/p_a(ij_com)

        ! derive interpolation weights
        call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
        call linint_Grow(ya_p, ya_l, ya_u, ya_grow, NY, iyl, iyr, varphi_y)

        ! get interpolated values of Qal and Qar
        if(varphi_a <= varphi_y)then
            Qal = varphi_a           *Q_al(ij_com, ial, iyl) + &
                  (varphi_y-varphi_a)*Q_al(ij_com, iar, iyl) + &
                  (1d0-varphi_y)     *Q_al(ij_com, iar, iyr)
            Qya = varphi_a           *Q_ya(ij_com, ial, iyl) + &
                  (varphi_y-varphi_a)*Q_ya(ij_com, iar, iyl) + &
                  (1d0-varphi_y)     *Q_ya(ij_com, iar, iyr)
        else
            Qal = varphi_y           *Q_al(ij_com, ial, iyl) + &
                  (varphi_a-varphi_y)*Q_al(ij_com, ial, iyr) + &
                  (1d0-varphi_a)     *Q_al(ij_com, iar, iyr)
            Qya = varphi_y           *Q_ya(ij_com, ial, iyl) + &
                  (varphi_a-varphi_y)*Q_ya(ij_com, ial, iyr) + &
                  (1d0-varphi_a)     *Q_ya(ij_com, iar, iyr)
        endif
        Qal = max(Qal, 1d-10)**(-1d0/gamma)
        Qya = max(Qya, 1d-10)**(-1d0/gamma)

        foc_annu = Qal/(Qya/p_a(ij_com)) - 1d0

    end function


    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = max(cons, 1d-10)**(-1d0/gamma)

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, ij, iy)

        implicit none
        integer, intent(in) :: ij, iy
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc, varphi, c_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij < JJ)then
            valuefunc = max(varphi*S(ij, ial, iy) + (1d0-varphi)*S(ij, iar, iy), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = c_help**egam/egam + beta*valuefunc

    end function


    ! get derivative of value function with respect to ya
    function dV_dya(a_plus, ij, iy)

        implicit none
        integer, intent(in) :: ij, iy
        real*8, intent(in) :: a_plus
        real*8 :: dV_dya, varphi_a
        integer :: ial, iar

        ! interpolate asset decision
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

        dV_dya = varphi_a*S_ya(ij, ial, iy) + (1d0-varphi_a)*S_ya(ij, iar, iy)
        dV_dya = beta**(-gamma)*dV_dya

    end function

end module
