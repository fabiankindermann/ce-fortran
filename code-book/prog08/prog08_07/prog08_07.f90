!##############################################################################
! PROGRAM EndogenousGrid
!
! ## Policy function iteration and the method of endogenous gridpoints
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
include "prog08_07m.f90"

program EndogenousGrid

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize a, c and value function
    call grid_Cons_Equi(a, 0d0, a0)
    c(:) = a(:)/2d0

    ! iterate until policy function converges
    do iter = 1, itermax

        ! set a = 0 manually
        a_endog(0) = 0d0
        c_endog(0) = 0d0

        ! calculate optimal decision for every gridpoint
        do ia = 1, NA

            ! calculate endogenous gridpoint and consumption
            c_endog(ia) = c(ia)*beta**(-gamma)
            a_endog(ia) = a(ia) + c_endog(ia)

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
        c_t(0) = spline_eval(a_t(0), coeff_c, 0d0, a0)
        do it = 1, TT
            a_t(it) = a_t(it-1) - c_t(it-1)
            c_t(it) = spline_eval(a_t(it), coeff_c, 0d0, a0)
        enddo
        call plot((/(dble(it),it=0,TT)/), c_t, legend='numerical')

        ! calculate the time path of consumption analytically
        a_t(0) = a0
        c_t(0) = a_t(0)*(1d0-beta**gamma)
        do it = 1, TT
            a_t(it) = a_t(it-1) - c_t(it-1)
            c_t(it) = a_t(it)*(1d0-beta**gamma)
        enddo
        call plot((/(dble(it),it=0,TT)/), c_t, legend='analytical')
        call execplot(xlabel='Time t', ylabel='Consumption c_t')


        ! plot numerical and analytical consumption
        call plot(a, c, legend='numerical')
        call plot(a, a*(1d0-beta**gamma), legend='analytical')
        call execplot(xlabel='Level of resources a', ylabel='Policy Function c(a)')

        ! calculate consumption function error
        err = 0d0
        do ia = 0, n_err
            a_err = a0*dble(ia)/dble(n_err)
            err_temp = abs(spline_eval(a_err, coeff_c, 0d0, a0) - a_err*(1d0-beta**gamma)) &
                                                /max(a_err*(1d0-beta**gamma), 1d-10)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(a, es15.7)')'Consumption function error:',err

        ! quit program
        stop

    end subroutine

end program
