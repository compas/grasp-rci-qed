program matrixelements_rcisettings
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use grasptest_testing
    use grasp_rciqed_rcisettings
    use grasp_rciqed_lib9290_init, only: lib9290_init_constants
    implicit none

    character(*), parameter :: jobname = "testjob"
    real(real64), parameter :: rtol = 1e-10_dp
    logical :: success = .true.
    type(rcisettings) :: settings

    call lib9290_init_constants

    call write_testsettings
    if(.not.read_settings_toml(jobname, settings)) then
        print '(a)', "read_settings_toml failed"
        success = .false.
    endif

    call test_isequal(success, "settings%Z", settings%Z, 77.0_dp, rtol)
    call test_isequal(success, "settings%atomic_mass_amu", settings%atomic_mass_amu, 100.0_dp, rtol)

    call test_isequal(success, "settings%breit_enabled", settings%breit_enabled, .true.)
    call test_isequal(success, "settings%nms_enabled", settings%nms_enabled, .false.)

    if(.not.success) then
        print *, "matrixelements_rcisettings: Tests failed."
        stop 1
    end if

contains

    subroutine write_testsettings
        use decide_C, only: LTRANS, LNMS, LSMS, LVP, LSE
        use def_C, only: Z, EMN, AUMAMU, FMTOAU
        use grid_C, only: RNT, H, N
        use npar_C, only: NPARM, PARM

        ! Set up the GRASP common blocks
        LTRANS=.true.
        LNMS=.false.
        LSMS=.true.
        LVP=.false.
        LSE=.true.

        Z = 77
        EMN = 100.0_dp / AUMAMU

        call write_settings_toml(jobname, 1)

    end subroutine write_testsettings

end program matrixelements_rcisettings
