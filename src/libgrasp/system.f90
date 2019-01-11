!> Routines and functions for interacting with the OS, filesystem etc.
module grasp_rciqed_system
    implicit none

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

end module grasp_rciqed_system
