!##############################################################################
! PROGRAM RamseyGov
!
! ## The Ramsey growth model with government activity
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
include "prog09_02m.f90"

program RamseyGov

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    c(:) = k(:)**alpha - delta*k(:)

    ! initialize government consumption
    gbar = gy*(alpha/(delta+1d0/beta-1d0))**(alpha/(1d0-alpha))

    ! iterate until policy function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call spline_interp(c, coeff_c)

        ! calculate decisions for every gridpoint
        do ik = 0, NK

            ! set starting value and communicate resource level
            x_in = c(ik)
            k_com = k(ik)

            ! find the optimal consumption level
            call fzero(x_in, foc, check)
            if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

            ! get optimal consumption and value function
            c_new(ik) = x_in

        enddo

        ! get convergence level
        con_lev = maxval(abs(c_new(:) - c(:))/max(abs(c(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        c = c_new
    enddo

    write(*,*)'No Convergence'

contains


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: TTplot = 800
        real*8 :: theta, UU, UU_0
        integer, parameter :: n_err = 10000
        real*8 :: err, c_err, err_temp

        ! end timer
        call toc()

        ! calculate the time path of consumption and capital numerically
        theta = 1d0/beta-1d0
        k_t(0) = (alpha/(theta+delta))**(1d0/(1d0-alpha))
        c_t(0) = spline_eval(k_t(0), coeff_c, k_l, k_u)
        y_t(0) = k_t(0)**alpha
        i_t(0) = y_t(0)-c_t(0)-gbar
        w_t(0) = (1d0-alpha)*k_t(0)**alpha
        r_t(0) = alpha*k_t(0)**(alpha-1d0)-delta
        b_t(0) = 0d0
        do it = 1, TT
            k_t(it) = (1d0-delta)*k_t(it-1) + i_t(it-1)
            c_t(it) = spline_eval(k_t(it), coeff_c, k_l, k_u)
            y_t(it) = k_t(it)**alpha
            i_t(it) = y_t(it)-c_t(it)-gbar
            w_t(it) = (1d0-alpha)*k_t(it)**alpha
            r_t(it) = alpha*k_t(it)**(alpha-1d0)-delta
            b_t(it) = by*k_t(it)**alpha
            tauw_t(it-1) = (gbar+(1d0+r_t(it-1))*b_t(it-1)-b_t(it)-tau_r*r_t(it-1)*(k_t(it-1)+b_t(it-1)))/w_t(it-1)
        enddo
        tauw_t(TT) = (gbar+r_t(TT)*b_t(TT)-tau_r*r_t(TT)*(k_t(TT)+b_t(TT)))/w_t(it-1)

        ! calculate the steady state utility level
        UU = c_t(TT)**egam/egam/(1d0-beta)

        ! calculate the utility level
        UU_0 = 0d0
        do it = 0, TT-1
            UU_0 = UU_0 + beta**(it)*c_t(it)**egam/egam
        enddo
        UU_0 = UU_0 + beta**TT*UU

        ! write output
        write(*,'(/a/)')'STEADY STATE VALUES AND WELFARE'
        write(*,'(a)')'   tau_w     K/Y     C/Y     I/Y     G/Y       C         U       U_0'
        write(*,'(6f8.2,2f10.3/)')(/tauw_t(TT), k_t(TT)/y_t(TT)/4d0, c_t(TT)/y_t(TT), &
            i_t(TT)/y_t(TT), gbar/y_t(TT)/)*100d0, c_t(TT), UU, UU_0

        call plot((/(dble(it),it=0,TTplot)/), k_t(0:TTplot))
        call execplot(xlabel='Time t', ylabel='Capital Stock k_t')

        call plot((/(dble(it),it=0,TTplot)/), c_t(0:TTplot))
        call execplot(xlabel='Time t', ylabel='Consumption c_t')

        call plot((/(dble(it),it=0,TTplot)/), r_t(0:TTplot))
        call execplot(xlabel='Time t', ylabel='Interest rate r_t')

        call plot((/(dble(it),it=0,TTplot)/), tauw_t(0:TTplot))
        call execplot(xlabel='Time t', ylabel='Wage tax tauw')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            c_err = spline_eval(k_com, coeff_c, k_l, k_u)
            err_temp = abs(foc(c_err)/c_err)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(/a, es15.7)')'Euler equation error:',err

        ! quit program
        stop

    end subroutine

end program
