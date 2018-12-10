module toml
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env
    implicit none

    type cpptoml_table
        type(c_ptr) :: p
    end type cpptoml_table

    interface

        function parse_file_(filename) bind(c, name="toml_parse_file")
            use, intrinsic :: iso_c_binding
            character(1, kind=c_char), intent(in) :: filename
            type(c_ptr) :: parse_file_
        end function parse_file_

        function contains_(table, key) bind(c, name="toml_contains")
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: table
            character(1, kind=c_char), intent(in) :: key
            logical(c_bool) :: contains_
        end function contains_

        function get_double_(table, key, value) bind(c, name="toml_get_double")
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: table
            character(1, kind=c_char), intent(in) :: key
            real(c_double), intent(out) :: value
            logical(c_bool) :: get_double_
        end function get_double_

        function get_int_(table, key, value) bind(c, name="toml_get_int")
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: table
            character(1, kind=c_char), intent(in) :: key
            integer(c_int), intent(out) :: value
            logical(c_bool) :: get_int_
        end function get_int_

        function get_bool_(table, key, value) bind(c, name="toml_get_bool")
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: table
            character(1, kind=c_char), intent(in) :: key
            logical(c_bool), intent(out) :: value
            logical(c_bool) :: get_bool_
        end function get_bool_

        function get_string_(table, key, value, length) bind(c, name="toml_get_string")
            use, intrinsic :: iso_c_binding
            type(c_ptr), intent(in), value :: table
            character(1, kind=c_char), intent(in) :: key(*)
            type(c_ptr), intent(out) :: value
            integer, intent(out) :: length
            logical(c_bool) :: get_string_
        end function get_string_

        subroutine string_delete_(string) bind(c, name="toml_string_delete")
            use, intrinsic :: iso_c_binding
            type(c_ptr), value, intent(in) :: string
        end subroutine string_delete_

    end interface

contains

    subroutine parse_file(filename, table)
        character(len=*), intent(in) :: filename
        type(cpptoml_table), intent(out) :: table

        table%p = parse_file_(cstr(filename))
    end subroutine parse_file

    function contains(table, key)
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        logical :: contains
        contains = contains_(table%p, cstr(key))
    end function contains

    function get_double(table, key, value)
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        real(real64), intent(out)  :: value
        logical :: get_double
        get_double = get_double_(table%p, cstr(key), value)
    end function get_double

    function get_integer(table, key, value)
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        integer, intent(out) :: value
        logical :: get_integer
        get_integer = get_int_(table%p, cstr(key), value)
    end function get_integer

    subroutine get_integer_default(table, key, default, value)
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        integer, intent(in) :: default
        integer, intent(out) :: value
        if(.not.get_integer(table, key, value)) then
            value = default
        endif
    end subroutine get_integer_default

    function get_logical(table, key, value)
        use iso_c_binding, only: c_bool
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        logical, intent(out) :: value
        logical :: get_logical

        logical(c_bool) :: value_

        ! We can't pass a logical to a logical(c_bool) argument directly, hence
        ! we'll need to have a temporary local variable for that.
        value_ = value
        get_logical = get_bool_(table%p, cstr(key), value_)
        value = value_
    end function get_logical

    subroutine get_logical_default(table, key, default, value)
        use iso_c_binding, only: c_bool
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        logical, intent(in) :: default
        logical, intent(out) :: value

        logical(c_bool) :: value_
        if(.not.get_logical(table, key, value)) then
            value = default
        endif
    end subroutine get_logical_default

    function get_string(table, key, value)
        use iso_c_binding, only: c_char, c_ptr
        type(cpptoml_table), intent(in) :: table
        character(len=*), intent(in) :: key
        character(:), allocatable, intent(out) :: value
        logical :: get_string

        integer :: length
        character(1, kind=c_char), pointer :: c_str(:)
        type(c_ptr) :: c_str_ptr
        integer :: i

        get_string = get_string_(table%p, cstr(key), c_str_ptr, length)
        if(get_string) then
            ! We are getting a raw C pointer from the C function and must type
            ! cast that into a c_char string here on the Fortran side.
            call c_f_pointer(c_str_ptr, c_str, [length])

            ! Here we allocate and assign, element by element, the string to the
            ! output string variable.
            allocate(character(length) :: value)
            do i = 1, length
                value(i:i) = c_str(i)
            enddo

            ! Call out to C again to clear the memory, as toml_get_string
            ! allocates memory for the `value` string.
            call string_delete_(c_str_ptr)
        endif
    end function get_string

    pure function cstr(str)
        character(len=*), intent(in) :: str
        character(len_trim(str)+1) :: cstr
        cstr = trim(str) // c_null_char
    end function cstr

end module toml
