!! Contains the global state of libgrasprci
module libgrasprci
    use, intrinsic :: iso_c_binding
    use grasp_rciqed_kinds, only: real64
    implicit none

    ! Error handling
    character(:), allocatable :: lasterror
    character(:), allocatable, target :: lasterror_cstr

    logical :: constants_initialized = .false.
    logical :: isodata_loaded = .false.
    logical :: csl_loaded = .false.
    logical :: orbitals_loaded = .false.
    logical :: mixing_loaded = .false.

    logical :: qedvp_initialized = .false.

    integer :: j2max

    type matrix1p
        real(real64), allocatable :: matrix(:,:)
    end type matrix1p

    type onescalar_cache
        integer :: ic, ir
        integer :: k, l
        real(real64), allocatable :: tshell(:)
    end type onescalar_cache

contains

    function grasp_libgraspci_matrix1p(matrix_cptr, i, j) bind(c)
        type(c_ptr), intent(in), value :: matrix_cptr
        integer(c_int), intent(in), value :: i, j
        real(c_double) :: grasp_libgraspci_matrix1p

        type(matrix1p), pointer :: matrix

        call c_f_pointer(matrix_cptr, matrix)
        grasp_libgraspci_matrix1p = matrix%matrix(i, j)
    end function grasp_libgraspci_matrix1p

    !> Accessor function for the `j2max` global variable in the `libgrasprci`
    !! module.
    function grasp_libgraspci_j2max() bind(c)
        use, intrinsic :: iso_c_binding
        implicit none
        integer(c_int) :: grasp_libgraspci_j2max
        grasp_libgraspci_j2max = j2max
    end

    function grasp_ci_onescalar(ic, ir) bind(c)
        use orb_C, only: NW
        use onescalar_I

        integer(c_int), intent(in), value :: ic, ir
        type(c_ptr) :: grasp_ci_onescalar

        type(onescalar_cache), pointer :: cache

        allocate(cache)
        allocate(cache%tshell(NW))
        cache%ic = ic
        cache%ir = ir
        call ONESCALAR(ic, ir, cache%k, cache%l, cache%tshell)
        grasp_ci_onescalar = c_loc(cache)
    end

    subroutine grasp_ci_onescalar_show(cache_cptr) bind(c)
        use orb_C, only: NW

        type(c_ptr), value :: cache_cptr

        type(onescalar_cache), pointer :: cache
        integer :: i

        call c_f_pointer(cache_cptr, cache)

        print '("ic=",i0,", ir=",i0)', cache%ic, cache%ir
        print '("k=",i0,", l=",i0)', cache%k, cache%l
        do i = 1,NW
            print *, cache%tshell(i)
        enddo
    end

    ! function grasp_ci_matrixelement_1p() bind(c)
    ! end function grasp_ci_matrixelement_1p

    function grasp_ci_matrixelement_1p_cached(cache_cptr, matrix_cptr) bind(c)
        use grasp_rciqed_cimatrixelements, only: qed_se

        type(c_ptr), value :: cache_cptr, matrix_cptr
        real(c_double) :: grasp_ci_matrixelement_1p_cached

        type(onescalar_cache), pointer :: cache
        type(matrix1p), pointer :: matrix

        call c_f_pointer(cache_cptr, cache)
        call c_f_pointer(matrix_cptr, matrix)

        grasp_ci_matrixelement_1p_cached = qed_se(matrix%matrix, cache%ic, cache%ir, cache%k, cache%l, cache%tshell)
    end function grasp_ci_matrixelement_1p_cached

    function grasp_ci_qedse_matrix(setype, matrix_ptr) bind(c)
        use orb_C, only: NW
        use grasp_rciqed_qed, only: qedse

        integer(c_int), intent(in), value :: setype
        type(c_ptr), intent(out) :: matrix_ptr
        integer(c_int) :: grasp_ci_qedse_matrix

        type(matrix1p), pointer :: sematrix

        if(setype < 1 .or. setype > 4) then
            lasterror = "grasp_ci_qedse_matrix: invalid setype specified"
            grasp_ci_qedse_matrix = 1
            return
        endif

        allocate(sematrix)
        allocate(sematrix%matrix(NW,NW))
        call qedse(setype, sematrix%matrix)
        matrix_ptr = c_loc(sematrix)

        grasp_ci_qedse_matrix = 0 ! no error
    end function grasp_ci_qedse_matrix

    function grasp_ci_init(isodata, csl, orbitals, mixing) bind(c)
        use, intrinsic :: iso_c_binding
        use grasp_c_interface
        use grasp_rciqed_lib9290_init, only: lib9290_init_constants, lib9290_init_rkco_gg
        use grasp_rciqed_lib9290_files, only: load_isodata, load_csl, load_orbitals, load_mixing
        use grasp_rciqed, only: init_rkintc
        implicit none

        character(kind=c_char), intent(in) :: isodata(1), csl(1), orbitals(1), mixing(1)
        integer(c_int) :: grasp_ci_init

        character(:), allocatable :: filename

        ! Initialize machine constants
        call lib9290_init_constants

        ! Load the input files
        call load_isodata(from_cstring(isodata))
        call load_csl(from_cstring(csl))
        call load_orbitals(from_cstring(orbitals))
        call load_mixing(from_cstring(mixing))

        ! Initialize angular coefficients (?) and DC integrals
        call lib9290_init_rkco_gg
        call init_rkintc(j2max)

        grasp_ci_init = 0
    end function grasp_ci_init

end module libgrasprci
