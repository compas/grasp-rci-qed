program matrixelements_hydrogenic
    use grasp_kinds, only: real64, dp
    use grasp_lib9290_init
    use g2k_lib92, only: lib92_init_csls, rcicommon_init
    use g2k_librci
    use grasptest_lib9290_setup
    use grasptest_lib9290_hydrogenic
    implicit none

    type(orbital_definition), parameter, dimension(*) :: orbitals = (/ &
        orbital_definition(1, -1), & ! 1s
        orbital_definition(2, -1), & ! 2s
        orbital_definition(2,  1), & ! 2p-
        orbital_definition(2, -2)  & ! 2p
    /)

    logical :: tests_passed = .true.
    integer :: n, k, l, tmpk, tmpl
    real(real64) :: vp_value
    character(:), allocatable :: testdata

    if(.not.getenv_allocating("RCIQED_TESTDATA", testdata)) then
        print '("ERROR: Unable to get RCIQED_TESTDATA environment variable")'
        error stop
    endif
    print '("RCIQED_TESTDATA=",a)', testdata

    call setup(74.0_dp, 183.91033628717801_dp) ! W-184
    call allocate_hydrogenic_orbitals(orbitals)

    call lib92_init_csls(testdata//"/oxygen.c")
    call lib9290_init_rkco_gg
    call rci_common_init ! TODO: There are two rci_common_init routines at the
                         ! moment. This one is more complete though.

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

    if(.not.tests_passed) then
        print *, "qed_flambaum_hydrogenic_test: Tests failed."
        stop 1
    end if

contains

    subroutine verify_dcb(ic, ir, reference)
        use grasp_cimatrixelements
        use grasp_kinds, only: real64, dp

        integer, intent(in) :: ir, ic
        real(real64), intent(in) :: reference

        real(real64), parameter :: tol = 1e-10_dp

        real(real64) :: hij = 0.0_dp

        hij = dirac_potential(ic, ir)
        hij = hij + coulomb(ic, ir)
        hij = hij + breit(ic, ir)

        print '(i3,i3,4es20.10)', ic, ir, hij, reference, reldiff(hij, reference), tol
        call check_tolerance("DCB", hij, reference, tol)
    end subroutine verify_dcb

    subroutine verify_dcbmsvp(ic, ir, reference)
        use grasp_cimatrixelements
        use grasp_kinds, only: real64, dp

        integer, intent(in) :: ir, ic
        real(real64), intent(in) :: reference

        real(real64), parameter :: tol = 1e-10_dp

        real(real64) :: hij = 0.0_dp

        hij = dirac_potential(ic, ir)
        hij = hij + coulomb(ic, ir)
        hij = hij + breit(ic, ir)
        hij = hij + qed_vp(ic, ir)
        hij = hij + nms(ic, ir)
        hij = hij + sms(ic, ir)

        print '(i3,i3,4es20.10)', ic, ir, hij, reference, reldiff(hij, reference), tol
        call check_tolerance("DCB", hij, reference, tol)
    end subroutine verify_dcbmsvp

    function reldiff(a, b)
        use grasp_kinds, only: real64
        real(real64), intent(in) :: a, b
        real(real64) :: reldiff
        reldiff = abs(a-b) / max(abs(a), abs(b))
    end function reldiff

    subroutine check_tolerance(which, a, b, relative_tolerance)
        use grasptest_testing, only: within_tolerance

        character(*), intent(in) :: which
        real(real64), intent(in) :: a, b, relative_tolerance
        real(real64) :: relative_difference

        if (.not.within_tolerance(a, b, relative_tolerance)) then
            relative_difference = abs(a-b) / max(abs(a), abs(b))
            print '("  Test failed: ",a," not within tolerance. Rel.diff: ",es12.5,", tol: ",es12.5)', &
                which, relative_difference, relative_tolerance
            tests_passed = .false.
        endif
    end subroutine check_tolerance

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
