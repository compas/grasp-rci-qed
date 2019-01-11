!> Various handy routines useful for writing unit tests.
module grasptest_testing
    implicit none

    interface test_isequal
        module procedure test_isequal_real64, test_isequal_logical
    end interface test_isequal

contains

    !> Calculate the relative difference of `a` and `b`.
    !!
    !! The relative difference is defined as:
    !!
    !! \f[
    !!   \frac{|a-b|}{\max(|a|, |b|)}
    !! \f]
    !!
    !! @param a,b Input values.
    !! @returns The relative difference of `a` and `b`.
    function reldiff(a, b)
        use grasp_rciqed_kinds, only: real64
        real(real64), intent(in) :: a, b
        real(real64) :: reldiff
        reldiff = abs(a-b) / max(abs(a), abs(b))
    end function reldiff

    !> Checks if the difference of `a` and `b` are within the tolerance relative
    !! to \f$\max(|a|,|b|)\f$.
    !!
    !! @param a,b Values to be checked.
    !! @param relative_tolerance Relative tolerance \f$\sigma\f$.
    !! @returns Whether \f$|a-b| / \max(|a|,|b|) < \sigma\f$.
    function within_tolerance(a, b, relative_tolerance)
        use grasp_rciqed_kinds, only: real64

        real(real64), intent(in) :: a, b, relative_tolerance
        logical :: within_tolerance
        real(real64) :: relative_difference

        relative_difference = abs(a-b) / max(abs(a), abs(b))
        if (relative_difference < relative_tolerance) then
            within_tolerance = .true.
        else
            within_tolerance = .false.
        endif
    end function within_tolerance

    !> Tests if two floating point values are the same withing the specified
    !! relative tolerance.
    !!
    !! Specifically, it tests that
    !!
    !! \f[
    !!     \frac{|a-b|}{\max(|a|, |b|)} < \sigma
    !! \f]
    !!
    !! where \f$\sigma\f$ is the specified relative tolerance.
    !!
    !! If the test fails, `test_passed` gets set to `.false.`, but the value is
    !! unchanged if the tests pass. This allows for the pattern where the test
    !! function is called multiple times before checking if any of the tests
    !! failed:
    !!
    !! ```
    !! call test_equality(success, "T1", a1, b1, rtol)
    !! call test_equality(success, "T2", a2, b2, rtol)
    !! if(.not.success) then
    !!     print *, "Test failures occurred."
    !! endif
    !! ```
    !!
    !! On failure, a message gets printed into the standard output that a test
    !! failed, which also includes the `which` string, allowing for the failed
    !! test to be identified easily.
    !!
    !! @param test_passed Gets set to `.false.` if the test fails, and is left
    !!   unchanged otherwise.
    !! @param which A string that identifies the test.
    !! @param a,b Values to be compared.
    !! @param relative_tolerance The relative tolerance \f$\sigma\f$.
    subroutine test_isequal_real64(test_passed, which, a, b, relative_tolerance)
        use grasp_rciqed_kinds, only: real64

        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        real(real64), intent(in) :: a, b, relative_tolerance
        real(real64) :: relative_difference

        if (.not.within_tolerance(a, b, relative_tolerance)) then
            relative_difference = abs(a-b) / max(abs(a), abs(b))
            print '("  Test failed: ",a," not within tolerance. Rel.diff: ",es12.5,", tol: ",es12.5)', &
                which, relative_difference, relative_tolerance
            test_passed = .false.
        endif
    end subroutine test_isequal_real64

    !> Tests if two logical values are equal.
    subroutine test_isequal_logical(test_passed, which, a, b)
        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        logical, intent(in) :: a, b

        if (.not.(a.eqv.b)) then
            print '("  Test failed: logical values differ for ",a," (a=",L1,", b=",L1,")")', &
                which, a, b
            test_passed = .false.
        endif
    end subroutine test_isequal_logical

end module grasptest_testing
