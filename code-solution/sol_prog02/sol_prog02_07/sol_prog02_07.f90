!##############################################################################
! PROGRAM TaxFunction
!
! ## Compute the optimal tax-rate with the help of polynomial interpolation
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
include "sol_prog02_07m.f90"

program TaxFunction

    use toolbox
    use globals

    implicit none

    ! declaration of variables
    integer :: i_max
    real*8 :: tau_opt, t_opt, tau_in, fret

    ! initialize data for interpolation
    tau_data = (/37d0, 42d0, 45d0/)
    t_data = (/198.875d0, 199.50d0, 196.875d0/)

    ! plot tax function with polynomial interpolation on equidistant grid
    call grid_Cons_Equi(tau_plot, tau_l, tau_u)
    t_plot = poly_interpol(tau_plot, tau_data, t_data)

    ! plot tax function
    call plot(tau_plot, t_plot)
    call execplot(xlabel='Tax rate',ylabel='Tax Revenue')

    ! find maximum tax-rate using the maxloc command
    i_max = maxloc(t_plot, 1) - 1
    tau_opt = tau_plot(i_max)
    t_opt = t_plot(i_max)

    ! print output
    write(*,'(a)') 'maxloc'
    write(*,'(a)') '------'
    write(*,'(a,f8.2, a)') 'Optimal tax rate: ', tau_opt, ' %'
    write(*,'(a, f8.2,a)') 'Max  tax revenue: ', t_opt, ' Mrd.'

    ! initialize initial guess for tau_in
    tau_in = 37.5d0

    ! call subroutine fminsearch
    call fminsearch(tau_in, fret, tau_l, tau_u, tax_func)

    ! print output
    write(*,'(/a)') 'fminsearch'
    write(*,'(a)') '----------'
    write(*,'(a,f8.2, a)') 'Optimal tax rate: ', tau_in, ' %'
    write(*,'(a, f8.2,a)') 'Max  tax revenue: ', -fret, ' Mrd.'

end program
