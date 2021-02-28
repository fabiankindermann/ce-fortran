!##############################################################################
! PROGRAM TaxFunction
!
! ## Evaluate german tax code for the income read from the console
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
program TaxFunction

    implicit none

    ! declaration of variables
    real*8 :: y, x, z, T, average, marginal

    ! read the income from console
    write(*,'(a,$)')'Type in the income: '
    read(*,*)y

    x = (y - 8130d0)/10000d0
    z = (y - 13469d0)/10000d0

    if(y < 8131d0)then
        T = 0d0
        average = 0d0
        marginal = 0d0
    elseif(y < 13470d0)then
        T = (933.70d0*x + 1400d0)*x
        average = T/y
        marginal = (1867.4d0*x + 1400d0)/10000d0
    elseif(y < 52882d0)then
        T = (228.74d0*z + 2397d0)*z + 1014d0
        average = T/y
        marginal = (457.48d0*z + 2397d0)/10000d0
    elseif(y < 250731d0)then
        T = 0.42d0*y - 8196d0
        average = T/y
        marginal = 0.42d0
    else
        T = 0.45d0*y - 15718d0
        average = T/y
        marginal = 0.45d0
    endif

    ! print output
    write(*,'(a, f14.6)')'tax burden:       ', T
    write(*,'(a, f14.6)')'average tax:      ', average*100d0
    write(*,'(a, f14.6)')'marginal tax rate:', marginal*100d0

end program
