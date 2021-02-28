!##############################################################################
! PROGRAM Plotgraphs3D
!
! ## How to plot graphs in 3D using the toolbox
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
program Plotgraphs3D

    use toolbox

    implicit none
    integer, parameter :: nplot1 = 80
    integer, parameter :: nplot2 = 50
    real*8 :: x(0:nplot1), y(0:nplot2), z(0:nplot1, 0:nplot2)
    real*8 :: a(0:nplot1), b(0:nplot1), c(0:nplot1)
    integer :: i1, i2

    ! Initialize x values
    do i1 = 0, nplot1
        x(i1) = -5d0 + 10d0/dble(nplot1)*dble(i1)
    enddo

    ! initialize y values
    do i2 = 0, nplot2
        y(i2) = -3d0 + 6d0/dble(nplot2)*dble(i2)
    enddo

    ! get z values
    do i1 = 0, nplot1
        do i2 = 0, nplot2
            z(i1, i2) = sin(x(i1))*cos(y(i2))
        enddo
    enddo

    ! call 3D plotting routine
    call plot3d(x, y, z)

    ! call 3D plotting routine with complete configuration
    call plot3d(x, y, z, color='black', linewidth=0.5d0, marker=2, &
               markersize=0.3d0, noline=.false., &
               xlim=(/-6d0, 6d0/), xticks=1.0d0, xlabel='x-Axis', &
               ylim=(/-3d0, 3d0/), yticks=0.5d0, ylabel='y-Axis', &
               zlim=(/-1d0, 1d0/), zticks=0.2d0, zlabel='z-Axis', &
               zlevel=0d0, surf=.true., surf_color=2, &
               transparent=.true., view=(/70d0, 50d0/), &
               title='sin(x)*cos(y)', filename='testplot', &
               filetype='eps', output='testdata')


    ! initialize spiral data
    do i1 = 0, nplot1
        c(i1) = 20d0*dble(i1)/dble(nplot1)
        a(i1) = sin(c(i1))
        b(i1) = cos(c(i1))
    enddo

    ! plot spiral
    call plot3d(a, b, c, linewidth=2d0)

end program
