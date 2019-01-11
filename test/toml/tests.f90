program tomltest
    use grasp_rciqed_kinds, only: real64, dp
    use grasptest_testing
    use grasp_rciqed_toml

    character(len=:), allocatable :: valuestr

    type(cpptoml_table) :: table
    logical :: tests_passed = .true.

    ! TODO: test that parsing didn't fail?
    call parse_file("test.toml", table)

    call test_contains("hamiltonian")
    call test_contains("nucleus.Z")
    call test_contains("nucleus.atomic_mass_amu")
    call test_contains("grid")
    call test_contains("grid.H")
    call test_contains("grid.RNT")
    call test_contains("grid.N")
    call test_contains_fail("hamiltoniaan")

    call test_get_double("nucleus.Z", 42.0_dp, 1e-15_dp)
    call test_get_double("nucleus.atomic_mass_amu", 88.234_dp, 1e-15_dp)
    call test_get_double("grid.H", 2e-3_dp, 1e-15_dp)
    call test_get_double("grid.RNT", 1e-3_dp, 1e-15_dp)
    call test_get_double_fail("grid.RNTz")

    call test_get_integer("grid.N", 500)
    call test_get_integer_default("grid.N", 200, 500)
    call test_get_integer_default("grid.Q", 42, 42)

    call test_get_logical("hamiltonian.breit", .true.)
    call test_get_logical("hamiltonian.nms", .false.)
    call test_get_logical_default("hamiltonian.breit", .false., .true.)
    call test_get_logical_default("hamiltonian.nms", .true., .false.)
    call test_get_logical_default("hamiltonian.nonexist", .true., .true.)
    call test_get_logical_default("hamiltonian.nonexist", .false., .false.)

    call test_get_string("nucleus.model", "point")

    ! Done. Report and exit non-zero exit code if need be.
    if(.not.tests_passed) then
        print '(a)', "TEST INFO: Some tests failed."
        stop 1
    endif
    print '(a)', "TEST INFO: All tests succeeded."

contains

    subroutine test_contains(key)
        character(len=*), intent(in) :: key
        if(.not.contains(table, key)) then
            print '("TEST ERROR: contains(table, ",a,") -> .false.")', key
            tests_passed = .false.
        endif
    end subroutine test_contains

    subroutine test_contains_fail(key)
        character(len=*), intent(in) :: key
        if(contains(table, key)) then
            print '("TEST ERROR: contains(table, ",a,") -> .true.")', key
            tests_passed = .false.
        endif
    end subroutine test_contains_fail

    subroutine test_get_double(key, reference, rtol)
        character(len=*), intent(in) :: key
        real(real64), intent(in) :: reference, rtol
        real(real64) :: value

        if(.not.get_double(table, key, value)) then
            print '("TEST ERROR: get_double(",a,") failed")', key
            tests_passed = .false.
        elseif(.not.within_tolerance(value, reference, rtol)) then
            print '("TEST ERROR: Bad value from get_double(",a,")")', key
            print '("   returned: ",e22.16)', value
            print '("   expected: ",e22.16)', reference
            print '("   rel.diff: ",e22.16)', reldiff(value, reference)
            print '("  tolerance: ",e22.16)', rtol
            tests_passed = .false.
        endif
    end subroutine test_get_double

    subroutine test_get_double_fail(key)
        character(len=*), intent(in) :: key
        real(real64) :: value
        if(get_double(table, key, value)) then
            print '("TEST ERROR: get_double(",a,") did not fail")', key
            tests_passed = .false.
        endif
    end subroutine test_get_double_fail

    subroutine test_get_integer(key, reference)
        character(len=*), intent(in) :: key
        integer, intent(in) :: reference
        integer :: value

        if(.not.get_integer(table, key, value)) then
            print '("TEST ERROR: get_integer(",a,") failed")', key
            tests_passed = .false.
        elseif(value /= reference) then
            print '("TEST ERROR: Bad value from get_integer(",a,")")', key
            print '("   returned: ",i0)', value
            print '("   expected: ",i0)', reference
            tests_passed = .false.
        endif
    end subroutine test_get_integer

    subroutine test_get_integer_default(key, default, reference)
        character(len=*), intent(in) :: key
        integer, intent(in) :: default, reference
        integer :: value

        call get_integer_default(table, key, default, value)
        if(value /= reference) then
            print '("TEST ERROR: Bad value from get_integer(",a,")")', key
            print '("   returned: ",i0)', value
            print '("   expected: ",i0)', reference
            tests_passed = .false.
        endif
    end subroutine test_get_integer_default

    subroutine test_get_logical(key, reference)
        character(len=*), intent(in) :: key
        logical, intent(in) :: reference
        logical :: value

        if(.not.get_logical(table, key, value)) then
            print '("TEST ERROR: get_logical(",a,") failed")', key
            tests_passed = .false.
        elseif(value.neqv.reference) then
            print '("TEST ERROR: Bad value from get_logical(",a,")")', key
            print '("   returned: ",L1)', value
            print '("   expected: ",L1)', reference
            tests_passed = .false.
        endif
    end subroutine test_get_logical

    subroutine test_get_logical_default(key, default, reference)
        character(len=*), intent(in) :: key
        logical, intent(in) :: default, reference
        logical :: value

        call get_logical_default(table, key, default, value)
        if(value.neqv.reference) then
            print '("TEST ERROR: Bad value from get_logical(",a,")")', key
            print '("   returned: ",L1)', value
            print '("   expected: ",L1)', reference
            tests_passed = .false.
        endif
    end subroutine test_get_logical_default

    subroutine test_get_string(key, reference)
        character(len=*), intent(in) :: key, reference
        character(:), allocatable :: value

        if(.not.get_string(table, key, value)) then
            print '("TEST ERROR: get_string(",a,") failed")', key
            tests_passed = .false.
        elseif(value /= reference) then
            print '("TEST ERROR: Bad value from get_string(",a,")")', key
            print '("   returned: ''",a,"''")', value
            print '("   expected: ''",a,"''")', reference
            tests_passed = .false.
        endif
    end subroutine test_get_string

end program tomltest
