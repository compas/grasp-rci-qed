program rci
    use parameter_def, only: NNNW
    use grasp_rciqed_kinds, only: real64, dp
    use grasptest_testing, only: test_isequal
    use getcid_I, only: enabled_orbitals
    !use grasp_rciqed_breit, only: breit_init
    implicit none
    integer :: i
    logical, dimension(1:0) :: empty_logical
    integer, dimension(1:0) :: empty_integer
    logical :: success = .true.
    logical :: enabledorbs(NNNW)

    ! Check that our testing routines actually work as intended
    call test_isequal(success, "isempty(empty_logical)", size(empty_logical), 0)
    call test_ieqal(success, "resize_false:1", resize_false(0, empty_logical), empty_logical)
    call test_ieqal(success, "resize_false:2", resize_false(3, (/ .true. /)), (/.true.,.false.,.false./))
    call test_ieqal(success, "resize_false:3", resize_false(4, (/ .false., .true. /)), (/.false.,.true.,.false.,.false./))

    ! Test enabled_orbital
    call enabled_orbitals(empty_integer, enabledorbs)
    call test_ieqal(success, "enabled_orbitals:1", enabledorbs, resize_false(NNNW, empty_logical))

    call enabled_orbitals((/ 1 /), enabledorbs)
    call test_ieqal(success, "enabled_orbitals:2", enabledorbs, resize_false(NNNW, (/ .true. /)))

    call enabled_orbitals((/ 2, 1 /), enabledorbs)
    call test_ieqal(success, "enabled_orbitals:3", enabledorbs, resize_false(NNNW, (/ .true., .true. /)))

    call enabled_orbitals((/ 5, 2, 3 /), enabledorbs)
    call test_ieqal(success, "enabled_orbitals:4", enabledorbs, resize_false(NNNW, (/ .false., .true., .true., .false., .true./)))

    if(.not.success) then
        print *, "breit: Tests failed."
        stop 1
    end if

contains

    !> Test Is EQual Array Logical
    subroutine test_ieqal(test_passed, which, a, b)
        use grasp_rciqed_kinds, only: real64

        logical, intent(inout) :: test_passed
        character(*), intent(in) :: which
        logical, intent(in) :: a(:), b(:)

        integer :: i
        logical :: arrays_equal = .true.

        if (size(a) /= size(b)) then
            print '("  Test failed: ",a," arrays are different length; size(a)=",i0,", size(b)=",i0)', &
                which, size(a), size(b)
            test_passed = .false.
        endif

        do i = 1, size(a)
            if(a(i) .neqv. b(i)) then
                arrays_equal = .false.
                exit
            endif
        enddo

        if (.not.arrays_equal) then
            print '("  Test failed: ",a," arrays are not equal at index ",i0)', which, i
            print *, "a=", a
            print *, "b=", b
            test_passed = .false.
        endif
    end subroutine test_ieqal

    function resize_false(nw, xs)
        integer, intent(in) :: nw
        logical, intent(in) :: xs(:)
        logical :: resize_false(nw)
        resize_false(1:nw) = .false.
        resize_false(1:size(xs)) = xs
    end

end program rci
