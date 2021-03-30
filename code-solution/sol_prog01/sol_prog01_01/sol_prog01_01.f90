!##############################################################################
! PROGRAM Calculation
!
! ## Arithmetic operations with "x" and "y" which are read from the console
!
! This code is published under the GNU General Public License v3
!                         (https://www.gnu.org/licenses/gpl-3.0.en.html)
!
! Authors: Hans Fehr, Maurice Hofmann and Fabian Kindermann
!          contact@ce-fortran.com
!
!##############################################################################
program Calculation

    implicit none

    ! declaration of variables
    integer :: ios
    real*8 :: x, y
    real*8 :: addit, diffit, prodit, quotit

    ! read the first number
    write(*,*)'Type in the 1. real number: '
    read(*,*, iostat = ios)x

    ! check whether input is readable
    if(ios /= 0)stop 'Error: number is not readable'

    ! read the second number
    write(*,*)'Type in the 2. real number: '
    read(*,*, iostat = ios)y

    ! check whether number is readable
    if(ios /= 0)stop 'Error: number is not readable'

    ! perform arithmetic operations
    addit = x + y
    diffit = x - y
    prodit = x*y
    quotit = x/y

    ! print output
    write(*,'(a, f14.6)')'Sum        =', addit
    write(*,'(a, f14.6)')'Difference =', diffit
    write(*,'(a, f14.6)')'Product    =', prodit
    write(*,'(a, f14.6)')'Quotient   =', quotit

end program
