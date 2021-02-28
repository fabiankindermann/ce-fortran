!##############################################################################
! PROGRAM oligopoly2
!
! ## Multi-dimensional rootfinding using the Gauss-Seidel iteration approach
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
program oligopoly2

    implicit none
    real*8 :: eta = 1.6d0
    real*8 :: c(2) = (/ 0.6d0, 0.8d0 /)
    real*8 :: qold(2), q(2), QQ, damp
    integer :: iter

    ! initialize q
    qold = 0.1d0
    damp = 0.7d0

    do iter = 1, 100

        QQ = sum(qold)
        q = 1d0/c*(QQ**(-1d0/eta)-1d0/eta*QQ**(-1d0/eta-1d0)*qold)
        q = damp*q + (1d0-damp)*qold

        ! write to screen
        write(*, '(a, i5,2f10.4)')'Iter: ', iter, q(1), q(2)

        ! check for convergence
        if( all(abs(q-qold) < 1d-6))then
            write(*,'(/a, 2f10.4)')'   Output', q(1), q(2)
            stop
        endif

        ! update q
        qold = q
    enddo

end program


