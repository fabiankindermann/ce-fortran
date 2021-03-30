!##############################################################################
! PROGRAM ConsumerSurplus
!
! ## Compute the the change in relative consumer surplus
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
include "sol_prog02_06m.f90"

program ConsumerSurplus

    use toolbox
    use globals

    implicit none

    ! declaration of variables
    integer :: i, j

    ! initialize consumption surplus
    CS = 0d0

    ! initialize lower bound for integral
    a = 0d0

    ! set up loop to compute the required intervals
    do i = 1, 2

        ! set p equal to the old market price
        if(i == 1)p = p_1

        ! set p equal to the new market price
        if(i == 2)p = p_2

        ! compute the upper bound for integral
        b = 2d0*p**(-0.5d0) - 1d0

        ! compute interval length
        h = (b-a)/dble(n)

        ! get equidistant nodes
        call grid_Cons_Equi(d, a, b)

        ! compute weights for trapezoid rule
        w_trapez(0) = h/2d0
        w_trapez(n) = h/2d0
        w_trapez(1:n-1) = h

        ! use trapezoid rule for integral
        f_trapez = (/(p_func(d(j)), j = 0, n)/)
        trapez = sum(w_trapez*f_trapez, 1)

        ! use simpson rule to compute integral
        simpson = p_func(a) + p_func(b)
        do j = 1, n-1
            simpson = simpson + 2d0*p_func(d(j))
        enddo
        do j = 0, n-1
            simpson = simpson + 4d0*p_func((d(j) + d(j+1))/2d0)
        enddo
        simpson = simpson*h/6d0

        ! use gauss-legendre to integral
        call legendre(a, b, d_gauss, w_gauss)
        f_gauss = (/(p_func(d_gauss(j)), j=0, n)/)
        gauss = sum(w_gauss*f_gauss, 1)

        ! compute consumer surplus
        CS(1, i) = trapez - b*p
        CS(2, i) = simpson - b*p
        CS(3, i) = gauss - b*p

    enddo

    ! print output
    write(*,'(a)') 'Consumer Surplus in %'
    write(*,'(a)') '---------------------'
    write(*,'(a,f12.6)') 'Trapezoid Rule:      ', (CS(1, 2)/CS(1, 1) - 1d0)*100d0
    write(*,'(a,f12.6)') 'Simpson Rule:        ', (CS(2, 2)/CS(2, 1) - 1d0)*100d0
    write(*,'(a,f12.6)') 'Gauss-Legendre:      ', (CS(3, 2)/CS(3, 1) - 1d0)*100d0

end program
