!##############################################################################
! PROGRAM EndogenousGridBorrowing
!
! ## Use endogenous gridpoints to solve a borrowing constraint problem
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
include "sol_prog08_09m.f90"

program EndogenousGridBorrowing

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! get lower and upper bound of grid
    a_l = max(-w/r, -a_borrow)
    a_u = w/r

    ! initialize a, c and policy function
    call grid_Cons_Equi(a, a_l, a_u)
    c(:) = (r*a(:)+w)/2d0
    ia_low = 0

    ! iterate until policy function converges
    do iter = 1, itermax

        ! calculate optimal decision for every gridpoint
        do ia = 0, NA

            ! calculate endogenous gridpoint and consumption at this point
            c_endog(ia) = c(ia)*(beta*(1d0+r))**(-gamma)
            a_endog(ia) = (a(ia) + c_endog(ia) - w)/(1d0+r)

        enddo

        ! add separate gridpoint or borrowing constraint
        a_endog(-1) = -a_borrow
        c_endog(-1) = -r*a_borrow + w

        ! choose gridpoints for interpolation
        do ia = NA, -1, -1

            ! find the first endogenous gridpoint with binding constraint
            if(a_endog(ia) <= -a_borrow)then
                ia_low = ia
                exit
            endif

        enddo

        ! stretch out to exogenous grid again
        do ia = 0, NA
            c_new(ia) = linint_Gen(a(ia), a_endog(ia_low:NA), c_endog(ia_low:NA), ia)
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

        implicit none
        integer :: it

        ! end timer
        call toc()

        ! interpolate policy function
        call spline_interp(c, coeff_c)

        ! calculate the time path of consumption with borr. constraint
        a_t(0) = a0
        c_t(0) = spline_eval(a_t(0), coeff_c, a_l, a_u)
        do it = 1, TT
            a_t(it) = (1d0+r)*a_t(it-1) + w - c_t(it-1)
            c_t(it) = spline_eval(a_t(it), coeff_c, a_l, a_u)
        enddo
        call plot((/(dble(it),it=0,TT)/), c_t, legend='with borr. constraint')

        ! calculate the time path of consumption without borr. constraint (analytically)
        a_t(0) = a0
        c_t(0) = (a_t(0)+w/r)*(1d0+r-(beta*(1d0+r))**gamma)
        do it = 1, TT
            a_t(it) = (1d0+r)*a_t(it-1) + w - c_t(it-1)
            c_t(it) = (a_t(it)+w/r)*(1d0+r-(beta*(1d0+r))**gamma)
        enddo

        call plot((/(dble(it),it=0,TT)/), c_t, legend='without borr. constraint')
        call execplot(xlabel='Time t', ylabel='Consumption c_t')

        ! plot consumption
        call plot(a, c, legend='with borr. constraint')
        call plot(a, (a+w/r)*(1d0+r-(beta*(1d0+r))**gamma), legend='without borr. constraint')
        call execplot(xlabel='Level of resources a', ylabel='Policy Function c(a)')

        ! quit program
        stop

    end subroutine

end program
