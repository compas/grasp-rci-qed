program matrixelements_hydrogenic
    use grasp_kinds, only: real64, dp
    use g2k_lib92, only: lib92_init_csls, lib92_init_rkco_gg, rcicommon_init
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
        error stop
    endif
    print '("RCI_TESTDATA=",a)', testdata

    call setup(74.0_dp)
    call allocate_hydrogenic_orbitals(orbitals)

    call lib92_init_csls(testdata//"/oxygen.c")
    call lib92_init_rkco_gg
    call rcicommon_init

    call calculate_hij(1, 1)
    call calculate_hij(1, 2)
    call calculate_hij(2, 2)
    call calculate_hij(3, 3)
    call calculate_hij(4, 4)
    call calculate_hij(4, 5)
    call calculate_hij(5, 5)

contains

    subroutine calculate_hij(ic, ir)
        use grasp_cimatrixelements

        integer, intent(in) :: ir, ic
        real(real64) :: hij_dp, hij_coulomb, hij_breit

        hij_dp = dirac_potential(ic, ir)
        hij_coulomb = coulomb(ic, ir)
        hij_breit = breit(ic, ir)

        print '("H(",i0,",",i0"):")', ic, ir
        print *, hij_dp, hij_coulomb, hij_breit
        print '(d30.16)', hij_dp + hij_coulomb + hij_breit
    end subroutine calculate_hij

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
