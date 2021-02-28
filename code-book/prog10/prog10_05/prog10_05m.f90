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

    use toolbox

    implicit none

    ! number of years the household retires
    integer, parameter :: JR = 45

    ! number of years the household lives
    integer, parameter :: JJ = 80

    ! number of white noise (zeta) shocks
    integer, parameter :: NW = 7

    ! number of transitory (epsilon) shocks
    integer, parameter :: NS = 7

    ! number of rate of return (vtheta) shocks
    integer, parameter :: NR = 7

    ! number of eps-vtheta shocks
    integer, parameter :: NSR = NS*NR

    ! number of points on the asset grid
    integer, parameter :: NA = 100

    ! number of points on the cash-on-hand-grid
    integer, parameter :: NX = 100

    ! number of points on the annuity grid
    integer, parameter :: NAR = NA

    ! household preference parameters
    real*8, parameter :: gamma = 0.10d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.96d0

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
    real*8, parameter :: a_u    = 500d0
    real*8, parameter :: a_grow = 0.07d0

    ! size of the annuity grid and loading factor
    real*8, parameter :: ar_l = a_l
    real*8, parameter :: ar_u = a_u
    real*8, parameter :: ar_grow = a_grow
    real*8, parameter :: xi = 0.00d0

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.07d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! lower and upper bound for X-grid
    real*8 :: X_l, X_u

    ! price of the annuity
    real*8 :: p_a

    ! lower and upper bound for the eta-grid
    real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w, eff(JJ), pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ), o_coh(JJ), yp_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_a(JJ), cv_o(JJ), cv_yp(JJ)

    ! individual variables
    real*8 :: X(0:NX), a(0:NA), ar(0:NAR), c(JJ, 0:NX, 0:NAR)
    real*8 :: a_plus(JJ, 0:NX, 0:NAR), omega_plus(JJ, 0:NA, 0:NAR), omegar_plus(0:NA)
    real*8 :: V(JJ, 0:NX, 0:NAR) = 0d0, Q(JJ, 0:NA, 0:NAR), S(0:NA)
    real*8 :: phi_X(JJ, 0:NX, 0:NAR), phi_a(JJ, 0:NA, 0:NAR)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA, 0:NAR), V_ar(JJ, 0:NX, 0:NAR) = 0d0, S_ar(0:NA) = 0d0
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
        real*8 :: foc_port, omega_p, R_port, X_p, c_p, varphi, dist
        integer :: ixl, ixr, iw, isr

        ! store portfolio share
        omega_p  = x_in
        foc_port = 0d0

        if(ij_com >= JR-1)then
            do isr = 1, NSR

                ! get future cash-on-hand and interpolate
                R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                X_p = R_port*a(ia_com) + pen(ij_com+1) + ar(ir_com)/p_a
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! calculate consumption and FOC
                c_p =       varphi*c(ij_com+1, ixl, ir_com) + &
                      (1d0-varphi)*c(ij_com+1, ixr, ir_com)
                c_p = max(c_p, 1d-10)
                foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*margu(c_p)
            enddo
        else
            do iw = 1, NW
                do isr = 1, NSR

                    ! get future cash-on-hand and interpolate
                    R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                    X_p = R_port*a(ia_com)/eps(isr) + w*eff(ij_com+1)*zeta(iw)
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_epsvtheta(isr)

                    ! calculate consumption and FOC
                    c_p =       varphi*c(ij_com+1, ixl, ir_com) + &
                          (1d0-varphi)*c(ij_com+1, ixr, ir_com)
                    c_p = max(c_p, 1d-10)
                    foc_port = foc_port + dist*(mu_r + vtheta(isr))*a(ia_com)*margu(eps(isr)*c_p)
                enddo
            enddo
        endif

    end function


    ! the first order condition regarding annuitization problem
    function foc_annu(x_in)

        implicit none
        real*8, intent(in) :: x_in
        real*8 :: foc_annu, omegar_p
        real*8 :: omega_p, R_port, X_p, al_p, ar_p, c_p, Var_p
        real*8 :: varphi_X, varphi_a, varphi_r, dist, Qal, Qar
        integer :: ixl, ixr, ial, iar, irl, irr, isr

        ! store retirement asset
        omegar_p = x_in

        ! determine future liquid wealth and retirement assets
        al_p = (1d0-omegar_p)*a(ia_com)
        ar_p = omegar_p*a(ia_com)

        ! derive interpolation weights
        call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
        call linint_Grow(ar_p, ar_l, ar_u, ar_grow, NAR, irl, irr, varphi_r)

        ! get optimal investment strategy
        if(varphi_a <= varphi_r)then
            omega_p = varphi_a           *omega_plus(JR-1, ial, irl) + &
                      (varphi_r-varphi_a)*omega_plus(JR-1, iar, irl) + &
                      (1d0-varphi_r)     *omega_plus(JR-1, iar, irr)
        else
            omega_p = varphi_r           *omega_plus(JR-1, ial, irl) + &
                      (varphi_a-varphi_r)*omega_plus(JR-1, ial, irr) + &
                      (1d0-varphi_a)     *omega_plus(JR-1, iar, irr)
        endif

        Qal = 0d0
        Qar = 0d0
        do isr = 1, NSR

            ! get future cash-on-hand and interpolate
            R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
            X_p = R_port*al_p + pen(JR) + ar_p/p_a
            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

            ! get distributional weight
            dist = dist_epsvtheta(isr)

            ! get future consumption
            if(varphi_X <= varphi_r)then
                c_p = varphi_X           *c(JR, ixl, irl) + &
                      (varphi_r-varphi_X)*c(JR, ixr, irl) + &
                      (1d0-varphi_r)     *c(JR, ixr, irr)
            else
                c_p = varphi_r           *c(JR, ixl, irl) + &
                      (varphi_X-varphi_r)*c(JR, ixl, irr) + &
                      (1d0-varphi_X)     *c(JR, ixr, irr)
            endif
            c_p = max(c_p, 1d-10)

            ! get future value of V_ar
            if(varphi_X <= varphi_r)then
                Var_p = varphi_X           *V_ar(JR, ixl, irl) + &
                        (varphi_r-varphi_X)*V_ar(JR, ixr, irl) + &
                        (1d0-varphi_r)     *V_ar(JR, ixr, irr)
            else
                Var_p = varphi_r           *V_ar(JR, ixl, irl) + &
                        (varphi_X-varphi_r)*V_ar(JR, ixl, irr) + &
                        (1d0-varphi_X)     *v_ar(JR, ixr, irr)
            endif
            Var_p = max(Var_p, 1d-10)**(-1d0/gamma)

            ! get Qal and Qar
            Qal = Qal + dist*R_port*margu(c_p)
            Qar = Qar + dist*(margu(c_p)/p_a + Var_p)
        enddo

        foc_annu = Qal/Qar - 1d0

    end function


    ! calculates marginal utility of consumption
    function margu(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu

        margu = max(cons, 1d-10)**(-1d0/gamma)

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, ij, ir)

        implicit none
        integer, intent(in) :: ij, ir
        real*8, intent(in) :: a_plus, cons
        real*8 :: valuefunc, varphi, c_help
        integer :: ial, iar

        ! check whether consumption or leisure are too small
        c_help = max(cons, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        if(ij == JR-1)then
            valuefunc = max(varphi*S(ial) + (1d0-varphi)*S(iar), 1d-10)**egam/egam
        elseif(ij < JJ)then
            valuefunc = max(varphi*Q(ij, ial, ir) + &
                    (1d0-varphi)*Q(ij, iar, ir), 1d-10)**egam/egam
        endif

        ! add todays part and discount
        valuefunc = c_help**egam/egam + beta*psi(ij+1)*valuefunc

    end function


    ! get derivative of value function with respect to ar
    function dV_dar(a_plus, ij, ir)

        implicit none
        integer, intent(in) :: ij, ir
        real*8, intent(in) :: a_plus
        real*8 :: dV_dar, varphi_a, varphi_X, omega_p, R_port, dist, X_p, c_p, Var_p
        integer :: ial, iar, ixl, ixr, isr

        ! interpolate asset decision
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi_a)

        ! get investment decision
        omega_p = varphi_a      *omega_plus(ij, ial, ir) + &
                  (1d0-varphi_a)*omega_plus(ij, iar, ir)

        dV_dar = 0d0

        do isr = 1, NSR

            ! get cash-on-hand
            R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
            X_p = R_port*a_plus + pen(ij+1) + ar(ir)/p_a

            ! derive weights
            call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

            ! get distributional weight
            dist = dist_epsvtheta(isr)

            ! derive V_y
            c_p = varphi_X      *c(ij+1, ixl, ir) + &
                  (1d0-varphi_X)*c(ij+1, ixr, ir)
            c_p = max(c_p, 1d-10)

            if(ij < JJ-1)then
                Var_p = varphi_X      *V_ar(ij+1, ixl, ir) + &
                        (1d0-varphi_X)*V_ar(ij+1, ixr, ir)
                Var_p = max(Var_p, 1d-10)**(-1d0/gamma)
            else
                Var_p = 0d0
            endif

            dV_dar = dV_dar + dist*(margu(c_p)/p_a + Var_p)
        enddo
        dV_dar = (beta*psi(ij+1)*dV_dar)**(-gamma)

    end function

end module
