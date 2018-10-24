!> Contains routines that can be used to set up the different global state
!! required by the lib9290 routines.
module grasptest_lib9290_setup
    implicit none
    ! implicit real*8 (A-H, O-Z)
    ! include 'parameters.def'

    ! COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS &
    !       /DEF1/EMN,IONCTY,NELEC,Z &
    !       /DEF2/C &
    !       /DEF4/ACCY,NSCF,NSIC,NSOLV &
    !       /DEF9/CVAC,PI &
    !       /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N &
    !       /NPAR/PARM(2),NPARM &
    !       /ORB2/NCF,NW,PNTRIQ &
    !       /ORB4/NP(NNNW),NAK(NNNW) &
    !       /ORB10/NH(NNNW) &
    !       /TATB/TA(NNN1),TB(NNN1),MTP

contains

    subroutine setup_constants
        use setcon_I
        use setmc_I

        call setcon
        call setmc ! machine-dependent parameters
    end subroutine setup_constants

    subroutine setup_grid
        use parameter_def
        use def_C
        use grid_C
        use tatb_C
        use setqic_I
        use radgrd_I

        H = 5.0D-2
        RNT = 2.0D-6
        HP = 0.0D0
        N = NNNP
        MTP = NNNP

        call setqic
        call radgrd

        ! Stolen from rci, based on the various examples of setting ACCY there
        ACCY = H**6
    end subroutine setup_grid

    subroutine setup_nucleus
        use def_C
        use npar_C
        use nucpot_I

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        print *, TENMAX,EXPMAX,EXPMIN,PRECIS
        print *, CVAC, PI

        Z = 18.0D0
        NPARM = 2
        PARM(1) = 6.8839456274865651D-005
        PARM(2) = 9.8905913700962641D-006

        call nucpot
    end subroutine setup_nucleus

    function kappa_to_string(kappa)
        integer, intent(in) :: kappa
        character, parameter, dimension(*) :: spec_notation = ['s','p','d','f','g','h','i','k','l','m','n']
        character(2) :: kappa_to_string

        if(kappa > 0) then
            kappa_to_string(1:1) = spec_notation(abs(kappa) + 1)
            kappa_to_string(2:2) = '-'
        else
            kappa_to_string(1:1) = spec_notation(abs(kappa))
            kappa_to_string(2:2) = ' '
        end if
    end function kappa_to_string

end module grasptest_lib9290_setup
