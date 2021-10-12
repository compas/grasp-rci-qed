!> Module to initialize the various parts of the `lib9290` global state.
!!
!! This primarily means populating the various `libmod` modules (derived from
!! the old COMMON blocks) with appropriate values.
module grasp_rciqed_lib9290_init
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    !> Initialize physical and unit-related constants and machine-dependent
    !! parameters.
    !!
    !! Calls `SETCON` and `SETMC`.
    subroutine lib9290_init_constants
        use setcon_I
        use setmc_I

        call setcon
        call setmc ! machine-dependent parameters
    end subroutine lib9290_init_constants

    !> Initializes the radial grid, suitable for the given nuclear charge
    !! `nuclear_z`.
    !!
    !! @param nuclear_z Charge of the nucleus.
    subroutine lib9290_init_grid(nuclear_z)
        ! Global state:
        use parameter_def, only: NNNP
        use def_C, only: ACCY
        use grid_C, only: RNT, H, HP, N
        use tatb_C, only: MTP
        ! lib9290 interfaces:
        use setqic_I
        use radgrd_I

        real(real64), intent(in) :: nuclear_z

        ! Taken from getcid.f90
        RNT = 2.0D-6 / nuclear_z
        H = 5.0D-2
        HP = 0.0D0
        N = NNNP
        ACCY = H**6

        MTP = NNNP ! TODO: is this actually necessary?

        call setqic
        call radgrd
    end subroutine lib9290_init_grid

    !> Sets up the nuclear parameters for a point nucleus with charge `nuclear_z`.
    !!
    !! @param nuclear_z Charge of the nucleus.
    subroutine lib9290_init_nucleus_pnc(nuclear_z)
        use def_C, only: CVAC, C, PI, TENMAX, EXPMAX, EXPMIN, PRECIS, Z
        use npar_C, only: NPARM, PARM
        use nucpot_I

        real(real64), intent(in) :: nuclear_z

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        ! Set the nucleus up as a point source for a specific Z
        Z = nuclear_z
        NPARM = 0

        call NUCPOT
    end subroutine lib9290_init_nucleus_pnc

    !> Sets up the nuclear parameters for a Fermi nucleus with charge `nuclear_z` and Fermi
    !! nuclear parameters `a` and `c`.
    !!
    !! @param nuclear_z Charge of the nucleus.
    !! @param a The a parameter of the Fermi charge distribution (in a.u.).
    !! @param c The c parameter of the Fermi charge distribution (in a.u.).
    subroutine lib9290_init_nucleus_fnc(nuclear_z, dist_a, dist_c)
        use def_C, only: CVAC, C, PI, TENMAX, EXPMAX, EXPMIN, PRECIS, Z
        use npar_C, only: NPARM, PARM
        use nucpot_I

        real(real64), intent(in) :: nuclear_z, dist_a, dist_c

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        ! Set the nucleus up as a point source for a specific Z
        Z = nuclear_z
        NPARM = 2
        PARM(1) = dist_c
        PARM(2) = dist_a

        call NUCPOT
    end subroutine lib9290_init_nucleus_fnc

    !> Sets up the nuclear mass, important for nuclear mass shifts.
    !!
    !! @param nuclear_mass Mass of the nucleus in atomic mass units.
    subroutine lib9290_init_nucleus_mass(nuclear_mass)
        use def_C, only: EMN, AUMAMU

        real(real64), intent(in) :: nuclear_mass

        EMN = nuclear_mass / AUMAMU
    end subroutine lib9290_init_nucleus_mass

    !> Initialize the global state necessary for the `rkco_gg` routine.
    !!
    !! Calls `ALCBUF` and `FACTT`.
    !!
    !! `rkco_gg` is defined in `librang90`.
    subroutine lib9290_init_rkco_gg
        use alcbuf_I
        use factt_I

        ! Set up the table of logarithms. From rci3mpi.f
        !
        ! This sets up the FACTS common block, which is used by the CORD routine
        ! that gets passed to RKCO_GG.
        call FACTT

        ! RKCO_GG uses the BUFFER common block in it's dependencies, but only
        ! re-allocates via ALCBUF. ALCBUF is actually decently documented in a
        ! sense -- it makes it relatively clear that you need to call ALCBUF(1)
        ! at some point. In rci_mpi it is done in SETHAM, before the main loop.
        ! So we'll do the allocation here.
        call ALCBUF(1)
    end subroutine lib9290_init_rkco_gg

end module grasp_rciqed_lib9290_init
