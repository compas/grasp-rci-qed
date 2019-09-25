!! Public functions of the libgrasp-rci shared library.

!>
!!
!! Note: calling this function again will invalidate the previously returned
!! pointer.
function grasp_error_string() bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    implicit none
    type(c_ptr) :: grasp_error_string
    lasterror_cstr = trim(lasterror)//c_null_char
    grasp_error_string = c_loc(lasterror_cstr)
end

subroutine grasp_ci_init_qedvp() bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci, only: qedvp_initialized
    use grasp_rciqed_qed, only: init_vacuum_polarization
    implicit none

    call init_vacuum_polarization
    qedvp_initialized = .true.
end subroutine grasp_ci_init_qedvp

!!
function grasp_initialize_constants() bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    !use grasp_c_interface
    use grasp_rciqed_lib9290_init, only: lib9290_init_constants
    implicit none

    integer(c_int) :: grasp_initialize_constants

    call lib9290_init_constants
    constants_initialized = .true.
    grasp_initialize_constants = 0

end function grasp_initialize_constants

!!
function grasp_load_isodata(filename_cstr) bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    use grasp_c_interface
    use grasp_rciqed_lib9290_files, only: load_isodata
    implicit none

    character(kind=c_char), intent(in) :: filename_cstr(1)
    integer(c_int) :: grasp_load_isodata

    character(:), allocatable :: filename

    if(.not.constants_initialized) then
        lasterror = "grasp_load_isodata: constants not initialized"
        grasp_load_isodata = 1
        return
    endif

    filename = from_cstring(filename_cstr)
    print *, filename
    call load_isodata(filename)

    grasp_load_isodata = 0
    return

end function grasp_load_isodata

!>
function grasp_load_csl(filename_cstr) bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    use grasp_c_interface
    use grasp_rciqed_lib9290_files, only: load_csl
    implicit none

    character(kind=c_char), intent(in) :: filename_cstr(1)
    integer(c_int) :: grasp_load_csl

    character(:), allocatable :: filename

    if(.not.constants_initialized) then
        lasterror = "grasp_load_isodata: constants not initialized"
        grasp_load_csl = 1
        return
    endif

    filename = from_cstring(filename_cstr)
    call load_csl(filename)

    grasp_load_csl = 0
    return

end function grasp_load_csl

!>
function grasp_load_orbitals(filename_cstr) bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    use grasp_c_interface
    use grasp_rciqed_lib9290_files, only: load_orbitals
    implicit none

    character(kind=c_char), intent(in) :: filename_cstr(1)
    integer(c_int) :: grasp_load_orbitals

    character(:), allocatable :: filename

    if(.not.constants_initialized) then
        lasterror = "grasp_load_isodata: constants not initialized"
        grasp_load_orbitals = 1
        return
    endif

    filename = from_cstring(filename_cstr)
    call load_orbitals(filename)

    grasp_load_orbitals = 0
    return

end function grasp_load_orbitals


!>
function grasp_load_mixing(filename_cstr) bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci
    use grasp_c_interface
    use grasp_rciqed_lib9290_files, only: load_mixing
    implicit none

    character(kind=c_char), intent(in) :: filename_cstr(1)
    integer(c_int) :: grasp_load_mixing

    character(:), allocatable :: filename

    if(.not.constants_initialized) then
        lasterror = "grasp_load_isodata: constants not initialized"
        grasp_load_mixing = 1
        return
    endif

    filename = from_cstring(filename_cstr)
    call load_mixing(filename)

    grasp_load_mixing = 0
    return

end function grasp_load_mixing

function grasp_init_rkco_gg() bind(c)
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_lib9290_init, only: lib9290_init_rkco_gg
    implicit none

    integer(c_int) :: grasp_init_rkco_gg

    call lib9290_init_rkco_gg
    grasp_init_rkco_gg = 0
end

function grasp_init_dcb() bind(c)
    use, intrinsic :: iso_c_binding
    use libgrasprci, only: j2max
    use grasp_c_interface
    !use grasp_rciqed_lib9290_init, only: lib9290_init_rkco_gg
    use grasp_rciqed, only: init_rkintc
    use grasp_rciqed_breit, only: init_breit
    implicit none

    integer(c_int) :: grasp_init_dcb

    !call lib9290_init_rkco_gg
    call init_rkintc(j2max)
    call init_breit(j2max)

    grasp_init_dcb = 0
end

function grasp_ncsfs() bind(c)
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_lib9290_csls, only: ncsfs_global
    implicit none
    integer(c_int) :: grasp_ncsfs
    grasp_ncsfs = ncsfs_global()
end

function grasp_ci_diracnuclear(i, j) bind(c)
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_cimatrixelements
    implicit none

    integer(c_int), intent(in), value :: i, j
    real(c_double) :: grasp_ci_diracnuclear

    grasp_ci_diracnuclear = dirac_potential(i, j)
end

function grasp_ci_coulomb(i, j) bind(c)
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_cimatrixelements
    implicit none

    integer(c_int), intent(in), value :: i, j
    real(c_double) :: grasp_ci_coulomb

    grasp_ci_coulomb = coulomb(i, j)
end

function grasp_orbital_grid(i) bind(c)
    use, intrinsic :: iso_c_binding
    use grid_C, only: N, R
    implicit none

    integer(c_int), intent(in), value :: i
    real(c_double) :: grasp_orbital_grid

    print *, i, R(1), R(2), R(i)

    grasp_orbital_grid = R(i)
end

function grasp_ci_qedvp(i, j) bind(c)
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_cimatrixelements
    implicit none

    integer(c_int), intent(in), value :: i, j
    real(c_double) :: grasp_ci_qedvp

    grasp_ci_qedvp = qed_vp(i, j)
end

function grasp_orbitals_nw() bind(c)
    use, intrinsic :: iso_c_binding
    use orb_C, only: NW
    implicit none
    integer(c_int) :: grasp_orbitals_nw
    grasp_orbitals_nw = NW
end

subroutine grasp_orbitals(nw_c, np_c, nak_c) bind(c)
    use, intrinsic :: iso_c_binding
    use orb_C, only: NW, NP, NAK
    implicit none

    integer(c_int), intent(out) :: nw_c
    integer(c_int), intent(out) :: np_c(NW), nak_c(NW)

    nw_c = NW
    np_c(1:NW) = NP(1:NW)
    nak_c(1:NW) = NAK(1:NW)
end
