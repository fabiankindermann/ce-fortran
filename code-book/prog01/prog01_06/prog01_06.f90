!##############################################################################
! PROGRAM Arithmetics
!
! ## Basic mathematical calculations
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
program Arithmetics

    ! declaration of variables
    implicit none
    real*8 :: a, b, c, d

    ! executable code
    a = 6d0
    b = 2.5d0
    c = exp(b) + a**2*sqrt(b)
    d = max(a,b)*sign(b, a)/mod(9d0,5d0) + abs(c)

    write(*, '(4f10.4)')a,b,c,d

end program
