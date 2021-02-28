!##############################################################################
! PROGRAM dynamic
!
! ## The dynamic solution to the cake eating problem
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
program dynamic

    use toolbox

    implicit none

    ! model parameters
    real*8, parameter :: gamma = 0.5d0
    real*8, parameter :: egam = 1d0-1d0/gamma
    real*8, parameter :: beta = 0.95d0
    real*8, parameter :: a0 = 100d0

    ! other variables
    integer :: it, ia

    integer, parameter :: TT = 200
    real*8 :: c_t(0:TT), a_t(0:TT)

    integer, parameter :: NA = 1000
    real*8 :: a(0:NA), c(0:NA), V(0:NA)

    ! calculate the time path of consumption
    a_t(0) = a0
    c_t(0) = a_t(0)*(1d0-beta**gamma)
    do it = 1, TT
        a_t(it) = a_t(it-1) - c_t(it-1)
        c_t(it) = a_t(it)*(1d0-beta**gamma)
    enddo

    call plot((/(dble(it),it=0,TT)/), c_t)
    call execplot(xlabel='Time t', ylabel='Consumption c_t')

    ! plot policy and value function
    do ia = 0, NA
        a(ia)  = 1d0 + dble(ia)/dble(NA)*a0
        c(ia) = a(ia)*(1d0-beta**gamma)
        V(ia) = (1d0-beta**gamma)**(-1d0/gamma)*a(ia)**egam/egam
    enddo

    ! plot policy function
    call plot(a, c)
    call execplot(xlabel='Level of resources a', ylabel='Policy Function c(a)')

    ! plot value function
    call plot(a, V)
    call execplot(xlabel='Level of Resources a', ylabel='Value Function V(a)')

end program
