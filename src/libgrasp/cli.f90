!> Contains routines related to command line input.
module grasp_rciqed_cli
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    !> Reads a numeric input from the user and checks that it is within the given
    !! range (`min <= n <= max`).
    !!
    !! It repeatedly asks the user until they give a valid value.
    !!
    !! Returns the value supplied by the user.
    !!
    function getoption(min, max)
        use iso_fortran_env, only: stdin => input_unit

        integer, intent(in) :: min, max
        integer :: getoption
        integer :: n

        ! TODO: Should it terminate after a certain number of tries? With an ERROR?
        do while (.true.)
            read(stdin,*) n
            if(min <= n .and. n <= max) then
                exit
            endif
            print '("Invalid choice (",i0,"). Must be ",i0,"-",i0,". Try again:")', n, min, max
        enddo
        getoption = n
    end function getoption

end module grasp_rciqed_cli
