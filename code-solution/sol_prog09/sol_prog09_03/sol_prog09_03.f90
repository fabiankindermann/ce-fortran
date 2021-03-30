!##############################################################################
! PROGRAM RamseyProd
!
! ## Ramsey model with two production sectors
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog09_03m.f90"

program RamseyProd

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize grid and policy function
    call grid_Cons_Equi(k, k_l, k_u)
    kc = k - (delta*k)**(1d0/alpha_i)

    ! iterate until policy function converges
    do iter = 1, itermax

        ! interpolate coefficients
        call spline_interp(kc, coeff_kc)

        ! calculate optimal decision for every gridpoint
        do ik = 0, NK

            ! initialize starting value and communicate resource level
            x_in = kc(ik)
            k_com = k(ik)

            ! find the optimal consumption level
            call fzero(x_in, foc, check)
            if(check)write(*,*)'ERROR IN ROOTFINDING PROCESS'

            ! get optimal consumption and value function
            kc_new(ik) = x_in

        enddo

        ! get convergence level
        con_lev = maxval(abs(kc_new(:) - kc(:))/max(abs(kc(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        kc = kc_new
    enddo

    write(*,*)'No Convergence'

contains


    ! For creating output plots
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: nplot = 1000
        real*8 :: kplot(0:nplot), kcplot(0:nplot)
        real*8 :: theta, kc_stat, ki_stat, k_stat, c_stat, i_stat
        integer, parameter :: n_err = 10000
        real*8 :: err, kc_err, err_temp

        ! end timer
        call toc()

        ! calculate stationary equilibrium values
        theta = 1d0/beta - 1d0
        ki_stat = (alpha_i/(theta + delta))**(1d0/(1d0 - alpha_i))
        k_stat = ki_stat**alpha_i/delta
        kc_stat = k_stat - ki_stat
        c_stat = kc_stat**alpha_c
        i_stat = ki_stat**alpha_i

        ! calculate the time path of consumption and capital numerically
        k_t(0) = k0
        kc_t(0) = spline_eval(k_t(0), coeff_kc, k_l, k_u)
        ki_t(0) = k_t(0) - kc_t(0)
        c_t(0) = kc_t(0)**alpha_c
        i_t(0) = ki_t(0)**alpha_i
        do it = 1, TT
            k_t(it) = (1d0-delta)*k_t(it-1) + ki_t(it-1)**alpha_i
            kc_t(it) = spline_eval(k_t(it), coeff_kc, k_l, k_u)
            ki_t(it) = k_t(it) - kc_t(it)
            c_t(it) = kc_t(it)**alpha_c
            i_t(it) = ki_t(it)**alpha_i
        enddo

        call plot((/(dble(it),it=0,TT)/), k_t)
        call plot((/0d0, dble(TT)/), (/k_stat, k_stat/))
        call execplot(xlabel='Time t', ylabel='Capital Stock k_t')

        call plot((/(dble(it),it=0,TT)/), kc_t)
        call plot((/0d0, dble(TT)/), (/kc_stat, kc_stat/))
        call execplot(xlabel='Time t', ylabel='Capital Stock kc_t')

        call plot((/(dble(it),it=0,TT)/), c_t, legend='Consumption')
        call plot((/0d0, dble(TT)/), (/c_stat, c_stat/))
        call plot((/(dble(it),it=0,TT)/), i_t, legend='Investment')
        call plot((/0d0, dble(TT)/), (/i_stat, i_stat/))
        call execplot(xlabel='Time t', ylabel='Consumption and Investment')

        ! calculate plot data for policy function
        do ik = 0, nplot
            kplot(ik) = k_l + (k_u-k_l)*dble(ik)/dble(nplot)
            kcplot(ik) = spline_eval(kplot(ik), coeff_kc, k_l, k_u)
        enddo

        ! plot policy function
        call plot(kplot, kcplot)
        call execplot(xlabel='Current level of capital k_t', ylabel='Policy Function kc(k_t)')

        ! calculate euler equation error
        err = 0d0
        do ik = 0, n_err
            k_com = k_l + (k_u-k_l)*dble(ik)/dble(n_err)
            kc_err = spline_eval(k_com, coeff_kc, k_l, k_u)
            err_temp = abs(foc(kc_err)/kc_err)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(a, es15.7)')'Euler equation error:',err

        ! quit program
        stop

    end subroutine

end program
