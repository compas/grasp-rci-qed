module grasptest_utilities
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    ! Borrowed from https://stackoverflow.com/a/31028207
    function itoa(i)
        character(len=:), allocatable :: itoa
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        itoa = trim(tmp)
    end

    function get_command_argument_allocating(n, value)
        logical :: get_command_argument_allocating
        integer, intent(in) :: n
        character(:), allocatable, intent(inout) :: value
        integer :: length, status
        character(1) :: test

        ! First, make an inquiry call to get_environment_variable to determine
        ! whether the variable exists and, if so, its length.
        call get_command_argument(n, test, length, status)
        if(status > 0) then
            ! From the GFortran manual:
            !
            ! > If the argument retrieval fails, STATUS is a positive number;
            !
            ! So a positive status will make this function fail (i.e. return
            ! .false.)
            get_command_argument_allocating = .false.
            return
        endif
        ! Will allocate or re-allocate value, unless it already has the correct length.
        if(allocated(value) .and. len(value) /= length) deallocate(value)
        if(.not.allocated(value)) allocate(character(length) :: value)
        ! Only call get_command_argument again if the argument is non-empty
        if(length > 0) then
            call get_command_argument(n, value, length, status)
        endif
        get_command_argument_allocating = .true.
    end

end module grasptest_utilities
