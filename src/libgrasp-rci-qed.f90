module libgrasp_rci_qed
    implicit none

contains

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

    subroutine libgrasp_rciqed_grid_c(n, rnt, h, hp, r_ptr, rp_ptr, rpor_ptr) bind(c)
        use, intrinsic :: iso_c_binding
        use grid_C, only: grid_c_n => N, grid_c_rnt => RNT, grid_c_h => H, grid_c_hp => HP, &
            grid_c_r => R, grid_c_rp => RP, grid_c_rpor => RPOR
        implicit none

        integer(c_int), intent(out) :: n
        real(c_double), intent(out) :: rnt, h, hp
        type(c_ptr), intent(out) :: r_ptr, rp_ptr, rpor_ptr

        real(c_double), pointer :: r(:), rp(:), rpor(:)

        n = grid_c_n
        rnt = grid_c_rnt
        h = grid_c_h
        hp = grid_c_hp

        allocate(r(n))
        r_ptr = c_loc(r)
        r(:) = grid_c_r(:n)

        allocate(rp(n))
        rp_ptr = c_loc(rp)
        rp(:) = grid_c_rp(:n)

        allocate(rpor(n))
        rpor_ptr = c_loc(rpor)
        rpor(:) = grid_c_rpor(:n)

    end subroutine libgrasp_rciqed_grid_c

    ! Initializes the the 9290 for a point nuclear charge with a given Z.
    ! This include initializing the physical constants, the grid and the nuclear
    ! potential itself.
    subroutine libgrasp_rciqed_9290_init_pnc(nuclear_z) bind(c)
        use, intrinsic :: iso_c_binding
        use grasp_rciqed_lib9290_init

        real(c_double), intent(in), value :: nuclear_z

        ! Calls SETCON and SETMC:
        call lib9290_init_constants
        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(nuclear_z)
        ! Set NPARM = 0 and calls NUCPOT:
        call lib9290_init_nucleus_pnc(nuclear_z)

    end subroutine libgrasp_rciqed_9290_init_pnc

    ! Initializes the the 9290 for a Fermi nucleus with a charge Z.
    ! This include initializing the physical constants, the grid and the nuclear
    ! potential itself.
    subroutine libgrasp_rciqed_9290_init_fnc(z, a, c) bind(c)
        use, intrinsic :: iso_c_binding
        use grasp_rciqed_lib9290_init

        real(c_double), intent(in), value :: z, a, c

        ! Calls SETCON and SETMC:
        call lib9290_init_constants
        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(z)
        ! Set NPARM = 1 and calls NUCPOT:
        call lib9290_init_nucleus_fnc(z, a, c)

    end subroutine libgrasp_rciqed_9290_init_fnc

    subroutine libgrasp_rciqed_vp_vacpol(array_n, zdist_ptr, vac2_ptr, vac4_ptr) bind(c)
        use grasp_rciqed_kinds, only: real64, dp
        use, intrinsic :: iso_c_binding
        ! Global state:
        use grid_C, only: R, RP
        use tatb_C, only: MTP, TB
        use ncdist_C, only: ZDIST
        ! Routines:
        use ncharg_I
        use vac2_I
        use vac4_I
        ! Arguments:
        integer(c_int), intent(out) :: array_n
        type(c_ptr), intent(out) :: zdist_ptr, vac2_ptr, vac4_ptr
        ! Local variables:
        real(c_double), pointer :: zdist_arr(:), vac2_arr(:), vac4_arr(:)

        ! Allocating the output arrays that go to the calling function:
        array_n = MTP
        allocate(zdist_arr(MTP))
        zdist_ptr = c_loc(zdist_arr)
        allocate(vac2_arr(MTP))
        vac2_ptr = c_loc(vac2_arr)
        allocate(vac4_arr(MTP))
        vac4_ptr = c_loc(vac4_arr)

        ! Initialize the nuclear charge distribution variable ZDIST. It shoud be all zeroes
        ! if it's PNC:
        call NCHARG
        ! Store the charge distribution:
        zdist_arr(:MTP) = ZDIST(:MTP)
        ! This "re-definition" of ZDIST is in VACPOL. As we call NCHARG every time anyway,
        ! it is fine to do this over and over again.
        ZDIST(:MTP) = ZDIST(:MTP)*R(:MTP)*RP(:MTP)

        ! VAC2, which should be calculating the Uehling potential, populates the TB array:
        call VAC2
        vac2_arr(:MTP) = TB(:MTP)

        ! VAC4, which should be calculating the Källen-Sabry potential, adds its potential
        ! to the TB array, on top of Uehling. So we need to subtract to get the original
        ! potential out.
        call VAC4
        vac4_arr(:MTP) = TB(:MTP) - vac2_arr(:MTP)
    end subroutine libgrasp_rciqed_vp_vacpol

    function libgrasp_rciqed_vp_funk(x, n) bind(c)
        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env, only: real64
        use funk_I

        real(c_double), intent(in), value :: x
        integer(c_int), intent(in), value :: n
        real(c_double) :: libgrasp_rciqed_vp_funk

        libgrasp_rciqed_vp_funk = funk(x, n)

    end function libgrasp_rciqed_vp_funk

    function libgrasp_rciqed_vp_funl(x, k) bind(c)
        use, intrinsic :: iso_c_binding
        use, intrinsic :: iso_fortran_env, only: real64
        use funl_I

        real(c_double), intent(in), value :: x
        integer(c_int), intent(in), value :: k
        real(c_double) :: libgrasp_rciqed_vp_funl

        libgrasp_rciqed_vp_funl = funl(x, k)

    end function libgrasp_rciqed_vp_funl

    ! This function tests calling subroutines that call subroutines in GRASP:
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

end module libgrasp_rci_qed
