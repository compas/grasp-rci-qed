!> Contains routines that can be used to set up the different global state
!! required by the lib9290 routines.
module grasptest_lib9290_setup
    implicit none

contains

    !> Call all the setup routines to set up parts of the lib9290 global state
    !! for a point nucleus with charge `nuclear_z`.
    subroutine setup(nuclear_z)
        use grasp_kinds, only: real64

        real(real64), intent(in) :: nuclear_z

        call setup_constants
        call setup_grid(nuclear_z)
        call setup_nucleus(nuclear_z)
    end subroutine setup

    subroutine setup_constants
        use setcon_I
        use setmc_I

        call setcon
        call setmc ! machine-dependent parameters
    end subroutine setup_constants

    subroutine setup_grid(nuclear_z)
        use grasp_kinds, only: real64
        use parameter_def
        use def_C
        use grid_C
        use tatb_C
        use setqic_I
        use radgrd_I

        real(real64), intent(in) :: nuclear_z

        ! Taken from getcid.f90
        RNT = 2.0D-6 / nuclear_z
        H = 5.0D-2
        HP = 0.0D0
        N = NNNP
        ACCY = H**6

        MTP = NNNP

        call setqic
        call radgrd
    end subroutine setup_grid

    !> Sets up the nuclear parameters for a point nucleus with charge `nuclear_z`.
    subroutine setup_nucleus(nuclear_z)
        use grasp_kinds, only: real64, dp
        use def_C, only: CVAC, C, PI, TENMAX, EXPMAX, EXPMIN, PRECIS, Z
        use npar_C, only: NPARM, PARM
        use nucpot_I

        real(real64), intent(in) :: nuclear_z

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        print *, TENMAX,EXPMAX,EXPMIN,PRECIS
        print *, CVAC, PI

        ! Set the nucleus up as a point source for a specific Z
        Z = nuclear_z
        NPARM = 0

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
