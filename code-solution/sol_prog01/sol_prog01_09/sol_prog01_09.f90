!##############################################################################
! PROGRAM CalcUtil
!
! ## Calculate utiltiy for different values of consumption and gamma
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
include "sol_prog01_09m.f90"

program CalcUtil

    use globals
    use toolbox

    implicit none

    ! declaration of variables
    integer :: ic, ig
    real*8 :: c_read, c(NC), u(NC)
    character(len = 100) :: label(NG)

    ! read in consumption value
    write(*,'(a, $)')'Type in a consumption level: '
    read(*,*)c_read

    ! print output
    write(*,*)
    do ig = 1, NG
        write(*,'(3(a,f8.2,2x))')'c = ', c_read, ' U = ', &
            utility(c_read, gamma_array(ig)), ' gamma = ', gamma_array(ig)
    enddo

    ! set up equidistant grid to plot utility to the screen
    do ic = 1, NC
        c(ic) = dble(ic)/dble(NC)
    enddo

    ! calculate and plot utility for different consumption levels
    do ig = 1, NG
        do ic = 1, NC
            u(ic) = utility(c(ic), gamma_array(ig))
        enddo

        ! print output
        write(label(ig),'(a,f6.4)')'gamma = ', gamma_array(ig)
        call plot(c, u, legend=label(ig))
        call execplot()
    enddo

end program
