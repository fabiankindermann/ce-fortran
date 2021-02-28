!##############################################################################
! PROGRAM creditriskmanagement_portfolio
!
! ## Credit risk management with a portfolio of bonds
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
program creditriskmanagement_portfolio

    implicit none
    integer, parameter :: KK = 8           ! number of ratings
    integer, parameter :: JJ = 5           ! number of industries
    integer, parameter :: BB = 16          ! number of bonds in portfolio
    real*8, parameter  :: sig_RR = 0.25d0  ! volatility of recovery rate
    real*8, parameter  :: r_f = 0.03d0     ! riskfree rate

    real*8 :: RM_data(KK-1, KK), cs_data(KK-1), RR_data(JJ)
    real*8 :: pi_D, NE, mu_RR, LGD, cs, cr, pi_k(KK)
    real*8 :: EL_d(BB+1), UL_d(BB+1), EL_m(BB+1), UL_m(BB+1)
    real*8 :: p_0, p_prime, sum_NE
    integer :: ij, ik, ik_p, it, TT, ib

    ! define a bond type
    type :: bond_type
        character(len=1) :: industry
        character(len=3) :: rating
        real*8 :: couponrate
        integer :: maturity
        real*8 :: exposure
    end type

    ! declare the bond variable
    type(bond_type) :: portfolio(BB)

    ! initialize bond data
    call init_data

    ! open file for output
    open(12, file = 'output.out')
    write(12,'(a)')'RISK ANALYSIS OF CORPORATE BOND PORTFOLIO'
    write(12, *)
    write(12,'(a)')'I  R        cr   T         p_0          NE        EL_D        UL_D        EL_M        UL_M'

    ! set up the portfolio
    portfolio(1)  = bond_type('F', 'BBB', 5.00d0,  4, 1000000d0)
    portfolio(2)  = bond_type('M', 'A'  , 4.63d0,  3, 3000000d0)
    portfolio(3)  = bond_type('M', 'BBB', 3.13d0,  6, 1000000d0)
    portfolio(4)  = bond_type('S', 'BBB', 2.88d0, 10, 1000000d0)
    portfolio(5)  = bond_type('T', 'BB' , 6.50d0,  2, 2000000d0)
    portfolio(6)  = bond_type('M', 'A'  , 5.00d0,  6, 5000000d0)
    portfolio(7)  = bond_type('M', 'A'  , 2.63d0,  5, 4000000d0)
    portfolio(8)  = bond_type('M', 'A'  , 1.88d0,  7,  500000d0)
    portfolio(9)  = bond_type('F', 'AA' , 6.50d0, 11, 8000000d0)
    portfolio(10) = bond_type('U', 'A'  , 6.38d0,  3, 4500000d0)
    portfolio(11) = bond_type('F', 'A'  , 4.00d0,  5, 1000000d0)
    portfolio(12) = bond_type('M', 'A'  , 2.00d0,  9,  600000d0)
    portfolio(13) = bond_type('M', 'A'  , 4.25d0,  5,  400000d0)
    portfolio(14) = bond_type('U', 'BB' , 5.13d0,  4, 2000000d0)
    portfolio(15) = bond_type('M', 'B'  , 4.38d0,  3, 1000000d0)
    portfolio(16) = bond_type('M', 'AA' , 5.63d0,  4,10000000d0)

    sum_NE = 0d0

    ! iterate over all bonds in portfolio
    do ib = 1, BB

        ! retrieve bond's properties
        cr    = portfolio(ib)%couponrate/100d0
        TT    = portfolio(ib)%maturity
        NE    = portfolio(ib)%exposure
        ij    = ij_code(portfolio(ib)%industry)
        ik    = ik_code(portfolio(ib)%rating)
        pi_D  = RM_data(ik, KK)
        pi_k  = RM_data(ik, :)
        cs    = cs_data(ik)
        mu_RR = RR_data(ij)

        ! compute price under current rating
        p_0 = 0d0
        do it = 1, TT-1
            p_0 = p_0 + cr/(1d0 + r_f + cs)**it
        enddo
        p_0 = p_0 + (1d0 + cr)/(1d0 + r_f + cs)**TT

        ! define loss given default
        LGD = p_0 - mu_RR
        sum_NE = sum_NE + NE

        ! default mode
        EL_d(ib) = NE*pi_D*LGD
        UL_d(ib) = NE*sqrt(pi_D*sig_RR**2 + LGD**2*pi_D*(1d0-pi_D))

        ! migration mode
        EL_m(ib) = 0d0
        UL_m(ib) = 0d0

        ! iterate over all potential credit events (except default)
        do ik_p = 1, KK-1

            cs = cs_data(ik_p)

            ! calculate price change if credit rating changes
            p_prime = 0d0
            do it = 1, TT-1
                p_prime = p_prime + cr/(1d0 + r_f + cs)**it
            enddo
            p_prime = p_prime + (1d0 + cr)/(1d0 + r_f + cs)**TT

            EL_m(ib) = EL_m(ib) + pi_k(ik_p)*(p_0 - p_prime)
            UL_m(ib) = UL_m(ib) + pi_k(ik_p)*(p_0 - p_prime)**2
        enddo

        ! add default state
        EL_m(ib) = EL_m(ib) + pi_D*LGD
        UL_m(ib) = UL_m(ib) + pi_D*(sig_RR**2 + LGD**2)

        ! determine EL and UL in migration mode
        EL_m(ib) = NE*EL_m(ib)
        UL_m(ib) = NE*sqrt(UL_m(ib) - (EL_m(ib)/NE)**2)

        write(12, '(a1,a5,f8.2,i4,f12.6,f12.0,4f12.2)')portfolio(ib)%industry, portfolio(ib)%rating, &
                        cr*100d0, TT, p_0, NE, EL_d(ib), UL_d(ib), EL_m(ib), UL_m(ib)
    enddo

    ! calculate expected and unexpected loss for total portfolio
    EL_d(BB+1) = sum(EL_d(1:BB))
    UL_d(BB+1) = sqrt(sum(UL_d(1:BB)**2))
    EL_m(BB+1) = sum(EL_m(1:BB))
    UL_m(BB+1) = sqrt(sum(UL_m(1:BB)**2))

    ! write total portfolio losses
    write(12,'(30x,f12.0,4f12.2)')sum_NE, EL_d(BB+1), UL_d(BB+1), EL_m(BB+1), UL_m(BB+1)
    write(12,'(44x,4(f10.2,a))')EL_d(BB+1)/sum_NE*10000d0, 'bp', UL_d(BB+1)/sum_NE*10000d0, 'bp', &
                                EL_m(BB+1)/sum_NE*10000d0, 'bp', UL_m(BB+1)/sum_NE*10000d0, 'bp'
    close(12)

contains


    ! initialize the probability parameters
    subroutine init_data()

        implicit none

        RM_data(1, :) = (/ 90.24d0,  8.99d0,  0.56d0,  0.05d0,  0.08d0,  0.03d0,  0.05d0,  0.00d0 /)
        RM_data(2, :) = (/  0.58d0, 90.03d0,  8.65d0,  0.56d0,  0.06d0,  0.08d0,  0.02d0,  0.02d0 /)
        RM_data(3, :) = (/  0.04d0,  2.00d0, 91.58d0,  5.71d0,  0.40d0,  0.17d0,  0.02d0,  0.08d0 /)
        RM_data(4, :) = (/  0.01d0,  0.13d0,  3.89d0, 90.69d0,  4.18d0,  0.68d0,  0.16d0,  0.26d0 /)
        RM_data(5, :) = (/  0.02d0,  0.04d0,  0.18d0,  5.81d0, 84.15d0,  7.97d0,  0.83d0,  1.00d0 /)
        RM_data(6, :) = (/  0.00d0,  0.05d0,  0.15d0,  0.25d0,  6.31d0, 83.16d0,  5.01d0,  5.07d0 /)
        RM_data(7, :) = (/  0.00d0,  0.00d0,  0.20d0,  0.30d0,  0.91d0, 15.96d0, 51.31d0, 31.32d0 /)
        RM_data = RM_data/100d0

        CS_data = (/0.00d0, 0.05d0, 1.00d0, 2.00d0, 5.00d0, 10.00d0, 20.00d0/)/100d0

        RR_data = (/0.437d0, 0.246d0, 0.575d0, 0.327d0, 0.493d0/)

    end subroutine init_data


    ! return rating code
    function ik_code(rating)

        implicit none
        character(len=3), intent(in) :: rating
        integer :: ik_code

        select case (rating)
            case('AAA')
                ik_code = 1
            case('AA')
                ik_code = 2
            case('A')
                ik_code = 3
            case('BBB')
                ik_code = 4
            case('BB')
                ik_code = 5
            case('B')
                ik_code = 6
            case('CCC')
                ik_code = 7
            case default
                stop 'Credit rating does not exists'
        end select

    end function


    ! return industry code
    function ij_code(industry)

        implicit none
        character(len=1), intent(in) :: industry
        integer :: ij_code

        select case (industry)
            case('M')
                ij_code = 1
            case('F')
                ij_code = 2
            case('U')
                ij_code = 3
            case('T')
                ij_code = 4
            case('S')
                ij_code = 5
            case default
                stop 'Industry does not exists'
        end select

    end function

end program

