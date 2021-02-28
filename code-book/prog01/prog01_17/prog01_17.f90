!##############################################################################
! PROGRAM Plotgraphs
!
! ## How to plot graphs with the plotting functions from the toolbox
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
program Plotgraphs

    use toolbox

    implicit none
    real*8 :: x(0:100), y(0:100)
    integer :: i1

    ! Initialize x values
    do i1 = 0, 100
        x(i1) = 1d0/100d0*dble(i1)
    enddo

    ! Calculate plot data
    y = x**2
    call plot(x, y)

    ! execute plot program
    call execplot()

    ! Calculate data for roots
    y = x**(1d0/2d0)

    ! you can specify a legend entry in the plot as follows
    call plot(x, y, legend='square root')

    ! the same for a cubic root
    y = x**(1d0/3d0)
    call plot(x, y, legend='cubic root')

    ! execute plot program and give the plot a title
    call execplot(title='Roots')

    ! plot has many more options that are specified here
    y = x**(1d0/2d0)
    call plot(x, y, color='green', linewidth=3d0, marker=2, &
            markersize=0.7d0,  noline=.false., legend='square root')

    y = x**(1d0/3d0)
    call plot(x, y, color='#5519D6', marker=5, markersize=1.2d0, &
              noline=.true., legend='cubic root')

    call execplot(xlim=(/0d0, 1.1d0/), xticks=0.1d0, &
           xlabel='x-Axis', ylim=(/0d0, 1.4d0/), yticks=0.2d0, &
           ylabel='y-Axis', title='Roots', legend='rs', &
           filename='testplot', filetype='eps', output='testdata')

end program
