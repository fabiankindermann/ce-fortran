!##############################################################################
! PROGRAM Sphere
!
! ## A program that uses a module
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
include "prog01_15m.f90"

program Sphere

    ! import variables and subroutines from module Volume
    use Volume

    implicit none

    write(*,'(f12.4)')vol(1d0)

end program
