!##############################################################################
! PROGRAM allinone
!
! ## The all-in-one solution to the cake eating problem
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
program allinone

    use toolbox

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: beta = 0.95d0
    real*8, parameter :: a0 = 100d0

    ! other variables
    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT)
    integer :: it

    ! calculate the time path of consumption
    do it = 0, TT
        c_t(it) = beta**(dble(it)*gamma)*(1d0 - beta**gamma)*a0
    enddo

    call plot((/(dble(it),it=0,TT)/), c_t)
    call execplot(xlabel='Time t', ylabel='Consumption c_t')

end program
