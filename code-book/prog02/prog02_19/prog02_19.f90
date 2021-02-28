!##############################################################################
! PROGRAM multi_interpol
!
! ## Multidimensional interpolation methods (linear and splines)
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
program multi_interpol

    use toolbox

    implicit none
    integer, parameter :: nx = 10, nz = 20, nerr = 1000
    real*8, parameter :: x_l = -1d0, x_r = 2d0
    real*8, parameter :: z_l = 0.5d0, z_r = 8d0
    real*8 :: xi(0:nx), zi(0:nz), yi(0:nx, 0:nz)
    real*8 :: xerr(0:nerr), zerr(0:nerr), yreal(0:nerr, 0:nerr)
    real*8 :: y_bi(0:nerr, 0:nerr), y_tr(0:nerr, 0:nerr), y_sp(0:nerr, 0:nerr)
    real*8 :: varphix, varphiz, coeff(nx+3, nz+3)
    integer :: ix, iz, ixl, ixr, izl, izr

    ! get nodes and data for interpolation
    call grid_Cons_Equi(xi, x_l, x_r)
    call grid_Cons_Equi(zi, z_l, z_r)
    do ix = 0, nx
        do iz = 0, nz
            yi(ix, iz) = xi(ix)**2 + sqrt(zi(iz))
        enddo
    enddo

    ! get nodes and data for calculation of errors
    call grid_Cons_Equi(xerr, x_l, x_r)
    call grid_Cons_Equi(zerr, z_l, z_r)
    do ix = 0, nerr
        do iz = 0, nerr
            yreal(ix, iz) = xerr(ix)**2 + sqrt(zerr(iz))
        enddo
    enddo

    ! piecewise linear interpolation
    do ix = 0, nerr
        do iz = 0, nerr
            call linint_Equi(xerr(ix), x_l, x_r, nx, ixl, ixr, varphix)
            call linint_Equi(zerr(iz), z_l, z_r, nz, izl, izr, varphiz)

            ! bilinear
            y_bi(ix, iz) = varphix*varphiz*yi(ixl, izl) &
                         + (1d0-varphix)*varphiz*yi(ixr, izl) &
                         + varphix*(1d0-varphiz)*yi(ixl, izr) &
                         + (1d0-varphix)*(1d0-varphiz)*yi(ixr, izr)

            ! on triangles
            if(varphix <= varphiz)then
                y_tr(ix, iz) = varphix*yi(ixl, izl) &
                             + (varphiz-varphix)*yi(ixr, izl) &
                             + (1d0-varphiz)*yi(ixr, izr)
            else
                y_tr(ix, iz) = varphiz*yi(ixl, izl) &
                             + (varphix-varphiz)*yi(ixl, izr) &
                             + (1d0-varphix)*yi(ixr, izr)
            endif
        enddo
    enddo

    ! cubic spline interpolation
    call spline_interp(yi, coeff)
    do ix = 0, nerr
        do iz = 0, nerr
            y_sp(ix, iz) = spline_eval((/xerr(ix), zerr(iz)/), &
                                coeff, (/x_l, z_l/), (/x_r, z_r/))
        enddo
    enddo

    ! write the maximum error to the screen
    write(*,'(a,f15.10)')'Maximum Error Bilinear:   ', &
                                maxval(abs(y_bi/yreal-1d0))*100d0
    write(*,'(a,f15.10)')'Maximum Error Triangular: ', &
                                maxval(abs(y_tr/yreal-1d0))*100d0
    write(*,'(a,f15.10)')'Maximum Error Spline:     ', &
                                maxval(abs(y_sp/yreal-1d0))*100d0

end program
