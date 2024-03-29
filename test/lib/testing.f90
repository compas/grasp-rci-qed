!> Various handy routines useful for writing unit tests.
module grasptest_testing
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

    interface test_isequal
        module procedure test_isequal_real64, test_isequal_logical, test_isequal_integer
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
    !! call test_isequal(success, "T1", a1, b1, rtol)
    !! call test_isequal(success, "T2", a2, b2, rtol)
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

        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        real(real64), intent(in) :: a, b, relative_tolerance
        real(real64) :: relative_difference

        if (.not.within_tolerance(a, b, relative_tolerance)) then
            relative_difference = abs(a-b) / max(abs(a), abs(b))
            print '("  Test failed: ",a," not within tolerance (value: ", &
                es12.5,", ref: ",es12.5,", rdiff: ",es12.5,", rtol: ",es12.5)', &
                which, a, b, relative_difference, relative_tolerance
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

    !> Tests if two integer values are equal.
    subroutine test_isequal_integer(test_passed, which, a, b)
        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        integer, intent(in) :: a, b

        if (a /= b) then
            print '("  Test failed: integer values differ for ",a," (a=",I0,", b=",I0,")")', &
                which, a, b
            test_passed = .false.
        endif
    end subroutine test_isequal_integer

    !> Tests if two floating point values are the same withing the specified
    !! absolute tolerance.
    !!
    !! Specifically, it tests that
    !!
    !! \f[
    !!     |a-b| < \sigma
    !! \f]
    !!
    !! where \f$\sigma\f$ is the specified absolute tolerance.
    !!
    !! If the test fails, `test_passed` gets set to `.false.`, but the value is
    !! unchanged if the tests pass. This allows for the pattern where the test
    !! function is called multiple times before checking if any of the tests
    !! failed:
    !!
    !! ```
    !! call test_isequal_atol(success, "T1", a1, b1, atol)
    !! call test_isequal_atol(success, "T2", a2, b2, atol)
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
    !! @param absolute_tolerance The absolute tolerance \f$\sigma\f$.
    subroutine test_isequal_atol(test_passed, which, a, b, absolute_tolerance)

        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        real(real64), intent(in) :: a, b, absolute_tolerance
        real(real64) :: absolute_difference

        absolute_difference = abs(b - a)
        if (absolute_difference > absolute_tolerance) then
            print '("  Test failed: ",a," not within abs. tolerance. Abs.diff: ",es12.5,", atol: ",es12.5)', &
                which, absolute_difference, absolute_tolerance
            test_passed = .false.
        endif
    end subroutine test_isequal_atol

end module grasptest_testing
