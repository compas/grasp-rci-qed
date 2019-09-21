!! Convenience functions for converting between the C and Fortran interfaces
module grasp_c_interface
    implicit none

contains

    function strlen(str)
        use, intrinsic :: iso_c_binding

        character(kind=c_char), intent(in) :: str(*)
        integer :: strlen

        strlen = 0
        do
            if(str(strlen+1) == c_null_char) then
                exit
            endif
            strlen = strlen + 1
        enddo
    end

    function from_cstring(str)
        use, intrinsic :: iso_c_binding

        character(kind=c_char), intent(in) :: str(*)
        integer :: length
        character(:), allocatable :: from_cstring

        length = strlen(str)
        allocate(character(length) :: from_cstring)
        from_cstring = transfer(str(1:length), from_cstring)
    end

end module grasp_c_interface
