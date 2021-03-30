!##############################################################################
! PROGRAM EndogenousGridINC
!
! ## Policy function iteration and the method of endogenous gridpoints
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog08_06m.f90"

program EndogenousGridINC

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! get lower and upper bound of grid
    a_l = -w/r
    a_u = w/r

    ! initialize a, c and value function
    call grid_Cons_Equi(a, a_l, a_u)
    c(:) = (r*a(:)+w)/2d0

    ! iterate until policy function converges
    do iter = 1, itermax

        ! set a = 0 manually
        a_endog(0) = a_l
        c_endog(0) = 0d0

        ! calculate optimal decision for every gridpoint
        do ia = 1, NA

            ! calculate endogenous gridpoint and consumption
            c_endog(ia) = c(ia)*(beta*(1d0+r))**(-gamma)
            a_endog(ia) = (a(ia) + c_endog(ia) - w)/(1d0+r)

        enddo

        ! stretch out to exogenous grid again
        do ia = 0, NA
            c_new(ia) = linint_Gen(a(ia), a_endog, c_endog, ia)
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

    ! For creating output plots.
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: n_err = 10000
        real*8 :: err, a_err, err_temp

        ! end timer
        call toc()

        ! interpolate policy function
        call spline_interp(c, coeff_c)

        ! calculate the time path of consumption numerically
        a_t(0) = a0
        c_t(0) = spline_eval(a_t(0), coeff_c, a_l, a_u)
        do it = 1, TT
            a_t(it) = (1d0+r)*a_t(it-1) + w - c_t(it-1)
            c_t(it) = spline_eval(a_t(it), coeff_c, a_l, a_u)
        enddo
        call plot((/(dble(it),it=0,TT)/), c_t, legend='numerical')

        ! calculate the time path of consumption analytically
        a_t(0) = a0
        c_t(0) = (a_t(0)+w/r)*(1d0+r-(beta*(1d0+r))**gamma)
        do it = 1, TT
            a_t(it) = (1d0+r)*a_t(it-1) + w - c_t(it-1)
            c_t(it) = (a_t(it)+w/r)*(1d0+r-(beta*(1d0+r))**gamma)
        enddo
        call plot((/(dble(it),it=0,TT)/), c_t, legend='analytical')
        call execplot(xlabel='Time t', ylabel='Consumption c_t')


        ! plot numerical and analytical consumption
        call plot(a, c, legend='numerical')
        call plot(a, (a+w/r)*(1d0+r-(beta*(1d0+r))**gamma), legend='analytical')
        call execplot(xlabel='Level of resources a', ylabel='Policy Function c(a)')

        ! calculate consumption function error
        err = a_l
        do ia = 0, n_err
            a_err = a_l + (a_u-a_l)*dble(ia)/dble(n_err)
            err_temp = abs(spline_eval(a_err, coeff_c, a_l, a_u) - &
                (a_err+w/r)*(1d0+r-(beta*(1d0+r))**gamma)) &
                /max((a_err+w/r)*(1d0+r-(beta*(1d0+r))**gamma), 1d-100)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(a, es15.7)')'Consumption function error:',err

        ! quit program
        stop

    end subroutine

end program
