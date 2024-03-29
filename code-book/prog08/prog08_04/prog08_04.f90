!##############################################################################
! PROGRAM GridSearch2
!
! ## Value function iteration and a grid search speed-up
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "prog08_04m.f90"

program GridSearch2

    use globals
    use toolbox

    implicit none

    ! start timer
    call tic()

    ! initialize a and value function
    call grid_Cons_Equi(a, 0d0, a0)
    V(:) = 0d0

    ! iterate until value function converges
    do iter = 1, itermax

        ! set a = 0 manually
        c(0) = 0d0
        ia_opt(0) = 0
        V_new(0) = V(1)-100d0

        ! calculate optimal decision for every gridpoint
        do ia = 1, NA

            ia_opt(ia) = max(ia-1, 0)
            cons = max(a(ia) - a(ia_opt(ia)), 1d-10)
            u_temp = cons**egam/egam + beta*V(ia_opt(ia))

            ! check whether there is an a_p that gives higher utility
            do ia_p = max(ia-2, 0), 0, -1

                ! calculate consumption and new utility
                cons = max(a(ia) - a(ia_p), 1d-10)
                u_new = cons**egam/egam + beta*V(ia_p)

                ! test which utility level is higher
                if(u_new > u_temp)then
                    ia_opt(ia) = ia_p
                    u_temp = u_new
                else
                    exit
                endif
            enddo

            ! get optimal consumption and value function
            c(ia) = a(ia)-a(ia_opt(ia))
            V_new(ia) = u_temp

        enddo

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
        real*8 :: err, err_temp

        ! end timer
        call toc()

        ! calculate the time path of consumption numerically
        ia_t(0) = NA
        c_t(0) = c(NA)
        do it = 1, TT
            ia_t(it) = ia_opt(ia_t(it-1))
            c_t(it) = c(ia_t(it))
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
        do ia = 0, NA
            err_temp = abs(c(ia) - a(ia)*(1d0-beta**gamma)) &
                                /max(a(ia)*(1d0-beta**gamma), 1d-10)
            if(err_temp > err)err = err_temp
        enddo
        write(*,'(a, es15.7)')'Consumption function error:',err

        ! quit program
        stop

    end subroutine

end program
