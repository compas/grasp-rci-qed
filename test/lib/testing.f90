!> Generic routines related useful in unit tests.
module grasptest_testing
    implicit none

contains

    !> Is the difference relative to `max(a,b)` within the tolerance value?
    function within_tolerance(a, b, relative_tolerance)
        use grasp_kinds, only: real64

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

end module grasptest_testing
