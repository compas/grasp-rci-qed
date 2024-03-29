program matrixelements_hydrogenic
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use grasp_rciqed_lib9290_init
    use grasp_rciqed_lib9290_files, only: load_csl
    use grasp_rciqed, only: init_rkintc
    use grasp_rciqed_breit, only: init_breit
    use grasp_rciqed_mass_shifts, only: init_mass_shifts
    use grasp_rciqed_qed_vp, only: qedvp_init
    use grasptest_lib9290_setup
    use grasptest_lib9290_hydrogenic
    implicit none

    type(orbital_definition), parameter, dimension(*) :: orbitals = (/ &
        orbital_definition(1, -1), & ! 1s
        orbital_definition(2, -1), & ! 2s
        orbital_definition(2,  1), & ! 2p-
        orbital_definition(2, -2)  & ! 2p
    /)

    real(real64), parameter :: rtol = 1e-10_dp

    logical :: success = .true.
    integer :: n, k, l, tmpk, tmpl, j2max
    real(real64) :: vp_value
    character(:), allocatable :: testdata

    if(.not.getenv_allocating("RCIQED_TESTDATA", testdata)) then
        print '("ERROR: Unable to get RCIQED_TESTDATA environment variable")'
        error stop
    endif
    print '("RCIQED_TESTDATA=",a)', testdata

    ! Setup for W-184
    call setup(74.0_dp, 183.91033628717801_dp)

    call allocate_hydrogenic_orbitals(orbitals)

    call load_csl(testdata//"/oxygen.c")
    call lib9290_init_rkco_gg

    ! RCI-specific initialization
    call init_rkintc(j2max)
    call init_breit(j2max)
    call qedvp_init
    call init_mass_shifts

    ! The following reference values were calculated using the non-MPI rci90
    ! program by adding in debug statements into setham_gg.f90. The values were
    ! taken from the EMT vector, with ELSTO added to all the diagonal values.
    !
    ! The physical system was a point-like Z=74 nucleus. The orbitals were
    ! generated using the exporthydrogenic tool and read from a file.
    !
    ! Note: the values here are calculated with orbitals that are re-generated
    ! with DCBSRW. Due to what are probably inaccuracies arising from
    ! interpolating the orbitals, the reference values do not exactly match the
    ! calculated values. However, it seems that they nicely fall into the 1e-10
    ! tolerance.
    !
    ! Also, any changes and errors in the routines populating hydrogenic orbitals
    ! may cause test failures here.
    print *
    print *, "Many-body DCB matrix elements for hydrogenic orbitals"
    print '(6(" "), 4a20)', &
        'H(ic,ir)', 'Reference', 'Diff (relative)', 'Tolerance'
    call verify_dcb(1, 1, -0.978959699717e4_dp)
    call verify_dcb(2, 2, -0.989063049881e4_dp)
    call verify_dcb(3, 3, -0.984231419490e4_dp)
    call verify_dcb(4, 4, -0.984161965244e4_dp)
    call verify_dcb(5, 5, -0.989314767535e4_dp)
    call verify_dcb(1, 2,  0.194514557433e1_dp)
    call verify_dcb(2, 1,  0.194514557433e1_dp)
    call verify_dcb(5, 4, -0.800396534942e0_dp)
    call verify_dcb(4, 5, -0.800396534942e0_dp)

    print *
    print *, "Many-body DCB + VP + NMS + SMS matrix elements for hydrogenic orbitals"
    print '(6(" "), 4a20)', &
        'H(ic,ir)', 'Reference', 'Diff (relative)', 'Tolerance'
    call verify_dcbmsvp(1, 1, -0.979221541410e4_dp)
    call verify_dcbmsvp(2, 2, -0.989328179926e4_dp)!
    call verify_dcbmsvp(3, 3, -0.984494905359e4_dp)!
    call verify_dcbmsvp(4, 4, -0.984425451113e4_dp)!
    call verify_dcbmsvp(5, 5, -0.989579897580e4_dp)!
    ! Note: the VP/NMS/SMS contributions to off-diagonal elements are zero.
    call verify_dcbmsvp(1, 2,  0.194514557433e1_dp)
    call verify_dcbmsvp(2, 1,  0.194514557433e1_dp)
    call verify_dcbmsvp(5, 4, -0.800396534942e0_dp)
    call verify_dcbmsvp(4, 5, -0.800396534942e0_dp)

    ! Mohr self-energy. Only diagonal elements can be estimated.
    print *
    print *, "Many-body SE (Mohr) contributions for hydrogenic orbitals"
    print '(6(" "), 4a20)', &
        'H(ic,ic)', 'Reference', 'Diff (relative)', 'Tolerance'
    call verify_se_mohr(1, 0.137580738458381e2_dp)
    call verify_se_mohr(2, 0.136535966904384e2_dp)
    call verify_se_mohr(3, 0.137058352681383e2_dp)
    call verify_se_mohr(4, 0.137058352681383e2_dp)
    call verify_se_mohr(5, 0.136535966904384e2_dp)

    if(.not.success) then
        print *, "matrixelements_hydrogenic: Tests failed."
        stop 1
    end if

contains

    subroutine verify_dcb(ic, ir, reference)
        use grasp_rciqed_cimatrixelements
        use grasptest_testing, only: reldiff, test_isequal

        integer, intent(in) :: ir, ic
        real(real64), intent(in) :: reference

        real(real64) :: hij = 0.0_dp

        hij = dirac_potential(ic, ir)
        hij = hij + coulomb(ic, ir)
        hij = hij + breit(ic, ir)

        print '(i3,i3,4es20.10)', ic, ir, hij, reference, reldiff(hij, reference), rtol
        call test_isequal(success, "DCB", hij, reference, rtol)
    end subroutine verify_dcb

    subroutine verify_dcbmsvp(ic, ir, reference)
        use grasp_rciqed_cimatrixelements
        use grasptest_testing, only: reldiff, test_isequal

        integer, intent(in) :: ir, ic
        real(real64), intent(in) :: reference

        real(real64) :: hij = 0.0_dp

        hij = dirac_potential(ic, ir)
        hij = hij + coulomb(ic, ir)
        hij = hij + breit(ic, ir)
        hij = hij + qed_vp(ic, ir)
        hij = hij + nms(ic, ir)
        hij = hij + sms(ic, ir)

        print '(i3,i3,4es20.10)', ic, ir, hij, reference, reldiff(hij, reference), rtol
        call test_isequal(success, "DCB", hij, reference, rtol)
    end subroutine verify_dcbmsvp

    subroutine verify_se_mohr(ic, reference)
        use grasp_rciqed_cimatrixelements
        use grasptest_testing, only: reldiff, test_isequal

        integer, intent(in) :: ic
        real(real64), intent(in) :: reference

        real(real64) :: hij = 0.0_dp

        hij = qed_se_mohr(ic)
        print '(2i3,4es20.10)', ic, ic, hij, reference, reldiff(hij, reference), rtol
        call test_isequal(success, "SEM", hij, reference, rtol)
    end subroutine verify_se_mohr

    !> Fetches an environment variable.
    !!
    !! If it was able to fetch the variable value, returns `.true.` and sets
    !! `value` to the value ([re]allocating if necessary). Returns `.false.` if
    !! the variable is not defined
    !!
    !! Under the hood it calls `get_environment_variable`, but it properly
    !! allocates or re-allocates the `value` to match the actual length of the
    !! environment variable.
    !!
    !! TODO: For this to be a proper library function, it should be implemented
    !! as an interface with additional methods to handle pointers and fixed-length
    !! strings.
    !!
    !! NOTE: It appears that this implementation does not work with GFortran 4.8
    !! and 4.9, probably due to a compiler bug (`value` in the caller context
    !! does not get set properly, even though it looks fine in
    !! `getenv_allocating`). It works with GFortran 5.5.0, but it is unknown if
    !! which 5.x version fixed the issue.
    function getenv_allocating(variable_name, value)
        logical :: getenv_allocating
        character(*), intent(in) :: variable_name
        character(:), allocatable, intent(inout) :: value
        integer :: length, status
        character(1) :: test

        ! First, make an inquiry call to get_environment_variable to determine
        ! whether the variable exists and, if so, its length.
        call get_environment_variable(variable_name, test, length, status)
        if(status > 0) then
            ! status of 1 or 2 from get_environment_variable means that the the
            ! variable is undefined, or that the system does not support
            ! environment variables at all, respectively. Both situations make
            ! getenv fail.
            getenv_allocating = .false.
            return
        endif
        ! Will allocate or re-allocate value, unless it already has the correct length.
        if(allocated(value) .and. len(value) /= length) deallocate(value)
        if(.not.allocated(value)) allocate(character(length) :: value)
        ! Can't pass a length 0 string to get_environment_variable.
        if(length > 0) then
            call get_environment_variable(variable_name, value, length, status)
        endif
        getenv_allocating = .true.
    end

end program matrixelements_hydrogenic
