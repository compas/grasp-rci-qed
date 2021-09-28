function libgrasp_rciqed_parameter_def(variable_name_cstr, variable_value) bind(c)
    use, intrinsic :: iso_c_binding
    use parameter_def
    implicit none

    character(kind=c_char), intent(in) :: variable_name_cstr(1)
    integer(c_int), intent(out) :: variable_value
    integer(c_int) :: libgrasp_rciqed_parameter_def

    character(:), allocatable :: variable_name

    variable_name = from_cstring(variable_name_cstr)
    if(variable_name == "KEYORB") then
        variable_value = KEYORB
    elseif(variable_name == "NNNP") then
        variable_value = NNNP
    elseif(variable_name == "NNN1") then
        variable_value = NNN1
    elseif(variable_name == "NNNW") then
        variable_value = NNNW
    elseif(variable_name == "NNNWM1") then
        variable_value = NNNWM1
    elseif(variable_name == "NNNWM2") then
        variable_value = NNNWM2
    else
        libgrasp_rciqed_parameter_def = 1
        return
    endif

    libgrasp_rciqed_parameter_def = 0 ! no error
    return

contains

    function from_cstring(str)
        use, intrinsic :: iso_c_binding

        character(kind=c_char), intent(in) :: str(*)
        integer :: length
        character(:), allocatable :: from_cstring

        length = strlen(str)
        allocate(character(length) :: from_cstring)
        from_cstring = transfer(str(1:length), from_cstring)
    end

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

end function libgrasp_rciqed_parameter_def


function libgrasp_rciqed_test(i1, i2, i3) bind(c)
    use, intrinsic :: iso_c_binding
    use itrig_I
    use setiso_I
    implicit none

    integer(c_int), intent(in), value :: i1, i2, i3
    integer(c_int) :: libgrasp_rciqed_test

    call setiso("isodata")

    libgrasp_rciqed_test = itrig(i1, i2, i3)

end function libgrasp_rciqed_test
