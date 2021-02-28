!##############################################################################
! PROGRAM creditriskmanagement
!
! ## Credit risk management with one bond
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
program creditriskmanagement

    implicit none
    integer, parameter :: KK = 8           ! number of ratings
    integer, parameter :: JJ = 5           ! number of industries
    real*8, parameter  :: sig_RR = 0.25d0  ! volatility of recovery rate
    real*8, parameter  :: r_f = 0.03d0     ! riskfree rate

    real*8 :: RM_data(KK-1, KK), cs_data(KK-1), RR_data(JJ)
    real*8 :: pi_D, NE, mu_RR, LGD, cs, cr, pi_k(KK)
    real*8 :: EL_d, UL_d, EL_m, UL_m
    real*8 :: p_0, p_prime
    integer :: ij, ik, ik_p, it, TT

    ! define a bond type
    type :: bond_type
        character(len=1) :: industry
        character(len=3) :: rating
        real*8 :: couponrate
        integer :: maturity
        real*8 :: exposure
    end type

    ! declare the bond variable
    type(bond_type) :: bond

    ! initialize bond data
    call init_data

    ! set up a bond
    bond = bond_type('F', 'BBB', 5.00d0,  4, 1000000d0)

    ! retrieve the bond's properties
    cr    = bond%couponrate/100d0
    TT    = bond%maturity
    NE    = bond%exposure
    ij    = ij_code(bond%industry)
    ik    = ik_code(bond%rating)
    pi_D  = RM_data(ik, KK)
    pi_k  = RM_data(ik, :)
    cs    = cs_data(ik)
    mu_RR = RR_data(ij)

    ! calculate price under current rating
    p_0 = 0d0
    do it = 1, TT-1
        p_0 = p_0 + cr/(1d0 + r_f + cs)**it
    enddo
    p_0 = p_0 + (1d0 + cr)/(1d0 + r_f + cs)**TT

    ! define loss given default
    LGD = p_0 - mu_RR

    ! default mode
    EL_d = NE*pi_D*LGD
    UL_d = NE*sqrt(pi_D*sig_RR**2 + LGD**2*pi_D*(1d0-pi_D))

    ! migration mode
    EL_m = 0d0
    UL_m = 0d0

    ! iterate over all potential credit events (except default)
    do ik_p = 1, KK-1

        cs = cs_data(ik_p)

        ! calculate price change if credit rating changes
        p_prime = 0d0
        do it = 1, TT-1
            p_prime = p_prime + cr/(1d0 + r_f + cs)**it
        enddo
        p_prime = p_prime + (1d0 + cr)/(1d0 + r_f + cs)**TT

        EL_m = EL_m + pi_k(ik_p)*(p_0 - p_prime)
        UL_m = UL_m + pi_k(ik_p)*(p_0 - p_prime)**2
    enddo

    ! add default state
    EL_m = EL_m + pi_D*LGD
    UL_m = UL_m + pi_D*(sig_RR**2 + LGD**2)

    ! determine EL and UL in migration model
    EL_m = NE*EL_m
    UL_m = NE*sqrt(UL_m - (EL_m/NE)**2)

    ! print results
    write(*,'(a)')'RISK ANALYSIS OF CORPORATE BOND'
    write(*, *)
    write(*,'(a)')'I  R        cr   T         p_0        EL_D        UL_D        EL_M        UL_M'
    write(*,'(a1,a5,f8.2,i4,f12.6,4f12.2)')bond%industry, bond%rating, &
                        cr*100d0, TT, p_0, EL_d, UL_d, EL_m, UL_m
    write(*,'(32x,4(f10.2,a))')EL_d/NE*10000d0, 'bp', UL_d/NE*10000d0, 'bp', &
                               EL_m/NE*10000d0, 'bp', UL_m/NE*10000d0, 'bp'

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
