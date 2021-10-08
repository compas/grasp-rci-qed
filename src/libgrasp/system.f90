!> Routines and functions for interacting with the OS, filesystem etc.
module grasp_rciqed_system
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

    !> Checks if `filename` exists or stops the program with `ERROR STOP`.
    !!
    !! The 3-argument version should be called by passing the `__FILE__` and
    !! `__LINE__` preprocessor variables as the second and third arguments:
    !!
    !! ```fortran
    !! call file_exists_or_stop(file_wfns, __FILE__, __LINE__)
    !! ```
    interface file_exists_or_stop
        module procedure file_exists_or_stop, file_exists_or_stop_where
    end interface file_exists_or_stop

contains

    !> Wraps the `inquire` statement to check if a file exists in a function.
    !!
    !! @param filename Name of the file to be checked.
    !! @return `.true.` if the file exists, `.false.` if not.
    function file_exists(filename)
        character(*), intent(in) :: filename
        logical :: file_exists
        inquire(file=filename, exist=file_exists)
    end function file_exists

    subroutine file_exists_or_stop(filename)
        character(*), intent(in) :: filename
        if(.not.file_exists(filename)) then
            print '("FATAL ERROR: File missing: ",a)', filename
            error stop
        endif
    end subroutine file_exists_or_stop

    subroutine file_exists_or_stop_where(filename, file, line)
        character(*), intent(in) :: filename, file
        integer, intent(in) :: line

        if(.not.file_exists(filename)) then
            print '("FATAL ERROR: File missing: ",a)', filename
            print '("Error raised at ",a,":",i0)', file, line
            error stop
        endif
    end subroutine file_exists_or_stop_where

end module grasp_rciqed_system
