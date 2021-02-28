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
    integer, parameter :: NAH = NA

    ! household preference parameters
    real*8, parameter :: gamma = 0.10d0
    real*8, parameter :: egam = 1d0 - 1d0/gamma
    real*8, parameter :: beta = 0.96d0
    real*8, parameter :: theta = 0.35d0
    real*8, parameter :: nu = 0d0   !10d0

    ! risk processes
    real*8, parameter :: sigma_zeta   = 0.0738d0
    real*8, parameter :: sigma_eps    = 0.0106d0
    real*8, parameter :: sigma_vtheta = 0.157d0**2d0
    real*8, parameter :: rho          = 0.00d0

    ! risk free rate and risk premium
    real*8, parameter :: r_f  = 0.02d0
    real*8, parameter :: mu_r = 0.04d0
    real*8, parameter :: load = 0.02d0

    ! size of the asset grid
    real*8, parameter :: a_l    = 0.0d0
    real*8, parameter :: a_u    = 450d0
    real*8, parameter :: a_grow = 0.07d0

    ! size of the housing grid
    real*8, parameter :: hmin = 0d0   !1d0
    real*8, parameter :: xi = 0.0d0   !0.7d0
    real*8, parameter :: delta_h = 0.01d0
    real*8, parameter :: ah_l = hmin
    real*8, parameter :: ah_u = a_u/(1d0-xi)
    real*8, parameter :: ah_grow = a_grow

    ! growth of the cash-on-hand grid
    real*8, parameter :: X_grow = 0.05d0

    ! pension fraction of last income
    real*8, parameter :: kappa = 0.50d0

    ! lower and upper bound for X-grid
    real*8 :: X_l, X_u

    ! lower and upper bound for the eta-grid
    real*8 :: eta_l(JJ), eta_u(JJ)

    ! discretized shocks
    real*8 :: dist_zeta(NW), zeta(NW)
    real*8 :: dist_epsvtheta(NSR), eps(NSR), vtheta(NSR)

    ! wages, transfer payments (old-age) and survival probabilities
    real*8 :: w, eff(JJ), pen(JJ), psi(JJ+1)

    ! cohort aggregate variables
    real*8 :: c_coh(JJ), y_coh(JJ), a_coh(JJ), ah_coh(JJ), o_coh(JJ)
    real*8 :: cv_c(JJ), cv_y(JJ), cv_a(JJ), cv_ah(JJ), cv_o(JJ)

    ! individual variables
    real*8 :: X(0:NX), a(0:NA), ah(0:NAH), c(JJ, 0:NX, 0:NAH)
    real*8 :: a_plus(JJ, 0:NX, 0:NAH), omega_plus(JJ, 0:NA, 0:NAH), omegah_plus(JJ, 0:NA)
    real*8 :: V(JJ, 0:NX, 0:NAH) = 0d0, Q(JJ, 0:NA, 0:NAH) = 0d0, S(JJ, 0:NA) = 0d0
    real*8 :: phi_X(JJ, 0:NX, 0:NAH), phi_a(JJ, 0:NA), phi_h(JJ, 0:NA, 0:NAH)

    ! numerical variables
    real*8 :: RHS(JJ, 0:NA)
    real*8 :: Q_al(1:JJ, 0:NA, 0:NAH), Q_ah(1:JJ, 0:NA, 0:NAH)
    integer :: ij_com, ix_com, ia_com, ih_com
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
        cons_com = X(ix_com) - a_plus - (1d0-xi)*hmin

        ! calculate linear interpolation for future part of first order condition
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        tomorrow = varphi*RHS(ij_com, ial) + (1d0-varphi)*RHS(ij_com, iar)

        ! calculate first order condition for consumption
        foc_cons = cons_com - tomorrow

    end function


    ! the first order condition regarding portfolio choice
    function foc_port(p)

        implicit none
        real*8, intent(in) :: p
        real*8 :: foc_port
        real*8 :: omega_p, R_port, X_p,  c_p, varphi_X, dist
        real*8 :: ah_p, varphi_h, beq
        integer :: ixl, ixr, ihl, ihr, iw, isr

        ! store portfolio share
        omega_p  = p

        foc_port = 0d0
        if(ij_com == JJ)then
            do isr = 1, NSR

                ! get future cash-on-hand and interpolate
                R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                if(a(ia_com)-xi*ah(ih_com) < 0d0)R_port = R_port + load

                ! get bequest
                beq = R_port*(a(ia_com)-xi*ah(ih_com)) + (1d0-delta_h)*ah(ih_com)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                foc_port = foc_port + dist*(mu_r + vtheta(isr))*(a(ia_com)-xi*ah(ih_com))* &
                    (1d0-psi(ij_com+1))*margu_b(beq)
            enddo
        elseif(ij_com >= JR-1)then
            do isr = 1, NSR

                ! get future cash-on-hand and interpolate
                R_port = 1d0 + r_f + omega_p*(mu_r + vtheta(isr))
                if(a(ia_com)-xi*ah(ih_com) < 0d0)R_port = R_port + load
                X_p = R_port*(a(ia_com)-xi*ah(ih_com)) + pen(ij_com+1) + (1d0-delta_h)*ah(ih_com)
                call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                ! get future house size
                ah_p = ah(ih_com)

                ! get bequest
                beq = R_port*(a(ia_com)-xi*ah(ih_com)) + (1d0-delta_h)*ah(ih_com)

                ! get distributional weight
                dist = dist_epsvtheta(isr)

                ! calculate consumption and FOC
                c_p = varphi_X      *c(ij_com+1, ixl, ih_com) + &
                      (1d0-varphi_X)*c(ij_com+1, ixr, ih_com)
                c_p = max(c_p, 1d-10)

                foc_port = foc_port + dist*(mu_r + vtheta(isr))*(a(ia_com)-xi*ah(ih_com))* &
                    (psi(ij_com+1)*margu_c(c_p) + (1d0-psi(ij_com+1))*margu_b(beq))
            enddo
        else
            do iw = 1, NW
                do isr = 1, NSR

                    ! get future cash-on-hand and interpolate
                    R_port = 1d0 + (r_f + omega_p*(mu_r + vtheta(isr)))
                    if(a(ia_com)-xi*ah(ih_com) < 0d0)R_port = R_port + load
                    X_p = R_port*(a(ia_com)-xi*ah(ih_com))/eps(isr) &
                          + w*eff(ij_com+1)*zeta(iw) + (1d0-delta_h)*ah(ih_com)/eps(isr)
                    call linint_Grow(X_p, X_l, X_u, X_grow, NX, ixl, ixr, varphi_X)

                    ! get future house size
                    ah_p = ah(ih_com)/eps(isr)
                    call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

                    ! get bequest
                    beq = R_port*(a(ia_com)-xi*ah(ih_com)) + (1d0-delta_h)*ah(ih_com)

                    ! get distributional weight
                    dist = dist_zeta(iw)*dist_epsvtheta(isr)

                    ! calculate consumption and FOC
                    if(varphi_X <= varphi_h)then
                        c_p = varphi_X           *c(ij_com+1, ixl, ihl) + &
                              (varphi_h-varphi_X)*c(ij_com+1, ixr, ihl) + &
                              (1d0-varphi_h)     *c(ij_com+1, ixr, ihr)
                    else
                        c_p = varphi_h           *c(ij_com+1, ixl, ihl) + &
                              (varphi_X-varphi_h)*c(ij_com+1, ixl, ihr) + &
                              (1d0-varphi_X)     *c(ij_com+1, ixr, ihr)
                    endif
                    c_p = max(c_p, 1d-10)

                    foc_port = foc_port + dist*(mu_r + vtheta(isr))*(a(ia_com)-xi*ah(ih_com))* &
                        (psi(ij_com+1)*margu_c(eps(isr)*c_p) + (1d0-psi(ij_com+1))*margu_b(beq))
                enddo
            enddo
        endif

    end function


    ! the first order condition regarding retirement savings problem
    function foc_real(p)

        implicit none
        real*8, intent(in) :: p
        real*8 :: foc_real
        real*8 :: omega_h, al_p, ah_p, Qal, Qar
        real*8 :: varphi_a, varphi_h
        integer :: ial, iar, ihl, ihr

        ! store share invested in retirement account
        omega_h  = p

        ! determine future liquid wealth and pre-return retirement wealth
        al_p = (1d0-omega_h)*a(ia_com)
        ah_p = hmin + omega_h*a(ia_com)/(1d0-xi)

        ! derive interpolation weights
        call linint_Grow(al_p, a_l, a_u, a_grow, NA, ial, iar, varphi_a)
        call linint_Grow(ah_p, ah_l, ah_u, ah_grow, NAH, ihl, ihr, varphi_h)

        ! get interpolated values of Qal and Qar
        if(varphi_a <= varphi_h)then
            Qal = varphi_a           *Q_al(ij_com, ial, ihl) + &
                  (varphi_h-varphi_a)*Q_al(ij_com, iar, ihl) + &
                  (1d0-varphi_h)     *Q_al(ij_com, iar, ihr)
            Qar = varphi_a           *Q_ah(ij_com, ial, ihl) + &
                  (varphi_h-varphi_a)*Q_ah(ij_com, iar, ihl) + &
                  (1d0-varphi_h)     *Q_ah(ij_com, iar, ihr)
        else
            Qal = varphi_h           *Q_al(ij_com, ial, ihl) + &
                  (varphi_a-varphi_h)*Q_al(ij_com, ial, ihr) + &
                  (1d0-varphi_a)     *Q_al(ij_com, iar, ihr)
            Qar = varphi_h           *Q_ah(ij_com, ial, ihl) + &
                  (varphi_a-varphi_h)*Q_ah(ij_com, ial, ihr) + &
                  (1d0-varphi_a)     *Q_ah(ij_com, iar, ihr)
        endif
        Qal = max(Qal, 1d-10)**(-1d0/gamma)
        Qar = max(Qar, 1d-10)**(-1d0/gamma)

        foc_real = Qal/(Qar/(1d0-xi)) - 1d0

    end function


    ! calculates marginal utility of consumption
    function margu_c(cons)

        implicit none
        real*8, intent(in) :: cons
        real*8 :: margu_c

        margu_c = theta*max(cons, 1d-10)**(-1d0/gamma)

    end function


    ! calculates marginal utility of housing
    function margu_h(house)

        implicit none
        real*8, intent(in) :: house
        real*8 :: margu_h

        margu_h = (1d0-theta)*max(house, 1d-10)**(-1d0/gamma)

    end function


    ! calculates marginal utility of bequest
    function margu_b(bequest)

        implicit none
        real*8, intent(in) :: bequest
        real*8 :: margu_b

        margu_b = nu*max(bequest, 1d-10)**(-1d0/gamma)

    end function


    ! calculates the value function
    function valuefunc(a_plus, cons, a_h, ij)

        implicit none
        integer, intent(in) :: ij
        real*8, intent(in) :: a_plus, cons, a_h
        real*8 :: valuefunc, varphi, c_help, ah_help
        integer :: ial, iar

        ! check whether consumption or house are too small
        c_help = max(cons, 1d-10)
        ah_help = max(a_h, 1d-10)

        ! get tomorrows utility
        call linint_Grow(a_plus, a_l, a_u, a_grow, NA, ial, iar, varphi)

        ! calculate tomorrow's part of the value function
        valuefunc = 0d0
        valuefunc = max(varphi*S(ij, ial) + (1d0-varphi)*S(ij, iar), 1d-10)**egam/egam

        ! add todays part and discount
        valuefunc = theta*c_help**egam/egam + (1d0-theta)*ah_help**egam/egam + beta*valuefunc

    end function

end module
