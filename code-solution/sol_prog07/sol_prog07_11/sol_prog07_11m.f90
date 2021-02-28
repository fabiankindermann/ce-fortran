!##############################################################################
! MODULE globals
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
module globals

    implicit none

    ! model parameters
    integer, parameter :: TT      = 36   ! 160
    integer, parameter :: JJ      = 6    !  55
    integer, parameter :: JR      = 5    ! 41
    real*8, parameter ::  gamma   = 0.5d0
    real*8, parameter ::  egam    = 1d0 - 1d0/gamma
    real*8, parameter ::  beta    = 0.9d0
    real*8, parameter ::  nubeq   = 3.0d0
    real*8, parameter ::  alpha   = 0.3d0
    real*8, parameter ::  delta   = 0.0d0
    real*8, parameter ::  tol     = 1d-6
    real*8, parameter ::  damp    = 0.3d0
    integer, parameter :: itermax = 250

    ! model variables
    real*8 :: w(0:TT), r(0:TT), wn(0:TT), Rn(0:TT), p(0:TT)
    real*8 :: tauw(0:TT), taur(0:TT), tauc(0:TT), taup(0:TT), pen(JJ,0:TT)
    real*8 :: by(0:TT), kappa(0:TT), n_p(0:TT), gy, tauk(0:TT)
    real*8 :: KK(0:TT), LL(0:TT), YY(0:TT), AA(0:TT), CC(0:TT), II(0:TT)
    real*8 :: BB(0:TT), GG(0:TT), BA(0:TT), BF(0:TT), TB(0:TT), BQ(0:TT)
    real*8 :: Tpen(0:TT), TXR(0:TT), h(JJ), m(JJ,0:TT)
    real*8 :: a(JJ+1,0:TT), c(JJ,0:TT), util(JJ,0:TT), v(-JJ+2:TT)
    real*8 :: psi(JJ+1,0:TT), beq(JJ,0:TT), omebeq(JJ), GAM(JJ,0:TT)
    integer :: tax(0:TT), eps(0:TT)
    integer :: it_com, ij_com
    logical :: lsra_on, smopec
    integer :: it

contains


    ! the first order condition
    function foc(x)

        implicit none
        real*8, intent(in) :: x(:)
        real*8 ::  foc(size(x, 1))
        real*8 :: vh, cons(ij_com:JJ), assets(ij_com:JJ+1)
        integer :: ij, itp, itm, ijp, it

        ! copy communication variable
        it = it_com
        ij = ij_com

        ! set up LSRA transfers
        vh = 0d0
        if(it == 1 .and. ij > 1)vh = v(2-ij)
        if(it >= 1 .and. ij == 1)vh = v(it)

        ! get assets
        if(ij == 1)then
            assets(1) = 0d0
        else
            assets(ij) = a(ij, it)
        endif

        ! get consumption and assets
        do ijp = ij, JJ

            itp = year(it, ij, ijp)
            cons(ijp) = x(ijp-ij+1)
            assets(ijp+1) = Rn(itp)*assets(ijp) + wn(itp)*h(ijp) + &
                beq(ijp,itp) + pen(ijp,itp) + vh - p(itp)*cons(ijp)
        enddo

        ! set up focs
        do ijp = ij, JJ

            itp = year(it, ij, ijp)
            itm = year(it, ij, ijp+1)

            if(ijp < JJ)then
                foc(ijp-ij+1) = beta*Rn(itm)*(psi(ijp+1,itm)*cons(ijp+1)**(-1d0/gamma)/p(itm) + &
                    (1d0-psi(ijp+1,itm))*nubeq*(Rn(itm)*assets(ijp+1))**(-1d0/gamma)) - &
                    cons(ijp)**(-1d0/gamma)/p(itp)
            else
                foc(ijp-ij+1) = beta*Rn(itm)*(nubeq*(Rn(itm)*assets(ijp+1))**(-1d0/gamma)) - &
                    cons(ijp)**(-1d0/gamma)/p(itp)
            endif
        enddo

    end function


    ! calculates year at which age ij agent is ijp
    function year(it, ij, ijp)

        implicit none
        integer, intent(in) :: it, ij, ijp
        integer :: year

        year = it + ijp - ij

        if(it == 0 .or. year <= 0)year = 0
        if(it == TT .or. year >= TT)year = TT

    end function

end module
