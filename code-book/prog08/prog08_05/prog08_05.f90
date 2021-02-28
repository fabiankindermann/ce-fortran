!##############################################################################
! PROGRAM Minimize
!
! ## Value function iteration and numerical minimization
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
include "prog08_05m.f90"

program Minimize

    use globals
    use toolbox

    implicit none
    integer :: ia, iter

    ! start timer
    call tic()

    ! initialize a, c and value function
    call grid_Cons_Equi(a, 0d0, a0)
    c(:) = a(:)/2d0
    V(:) = 0d0
    coeff_v(:) = 0d0

    ! iterate until value function converges
    do iter = 1, itermax

        ! set a = 0 manually
        c(0) = 0d0
        V_new(0) = -1d10

        ! calculate optimal decision for every gridpoint
        do ia = 1, NA

            ! initialize starting value and communicate resources
            x_in = a(ia) - c(ia)
            ia_com = ia

            call fminsearch(x_in, fret, 0d0, a(ia), utility)

            ! get optimal consumption and value function
            c(ia) = a(ia) - x_in
            V_new(ia) = -fret

        enddo

        ! interpolate coefficients
        call spline_interp((egam*V_new)**(1d0/egam), coeff_V)

        ! get convergence level
        con_lev = maxval(abs(V_new(:) - V(:))/max(abs(V(:)), 1d-10))
        write(*,'(i5,2x,f20.7)')iter, con_lev

        ! check for convergence
        if(con_lev < sig)then
            call output()
        endif

        V = V_new
    enddo

    write(*,*)'No Convergence'

contains


    ! For creating output plots.
    subroutine output()

        use toolbox

        implicit none
        integer, parameter :: n_err = 10000
        integer :: it
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

        ! plot numerical and analytical value function
        call plot(a(10:NA), V(10:NA), legend='numerical')
        call plot(a(10:NA), (1d0-beta**gamma)**(-1d0/gamma)*a(10:NA) &
                                                **egam/egam, legend='analytical')
        call execplot(xlabel='Level of Resources a', ylabel='Value Function V(a)')

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
