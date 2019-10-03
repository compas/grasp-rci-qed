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

    !> If set to `.true.`, will enable various debug output to the standard
    !! output.
    logical :: debug = .false.

    !> Stores the matrix elements of a scalar single-particle operator.
    type matrix1pscalar
        !> @var n
        !! The number \f$n\f$ of single particle orbitals forming the basis in
        !! which the operator is represented.
        integer :: n
        !> @var matrix
        !! An \f$n \times n\f$ matrix of the matrix elements.
        real(real64), allocatable :: matrix(:,:)
    end type matrix1pscalar

    !> Caches the values from the `ONESCALAR` routine, forming the energy
    !! expression of a 1-particle scalar operator between the CSFs `ic` and
    !! `ir`.
    type onescalar_cache
        !> @var ir
        !! The row index in the list of CSFs.
        !! @var ic
        !! The column index in the list of CSFs.
        integer :: ir, ic
        !> @var k
        !! The output `k` value from `ONESCALAR`.
        !! @var l
        !! The output `l` value from `ONESCALAR`.
        integer :: k, l
        !> @var tshell
        !! The output `tshell` array from `ONESCALAR`.
        real(real64), allocatable :: tshell(:)
    end type onescalar_cache

contains

    !===========================================================================
    ! The C API for managing the library global state
    !---------------------------------------------------------------------------

    !> Initializes the various global constants.
    subroutine libgrasprci_initialize_constants() bind(c)
        use grasp_rciqed_lib9290_init, only: lib9290_init_constants
        call lib9290_init_constants
        constants_initialized = .true.
    end subroutine libgrasprci_initialize_constants

    !> Accessor function for the `j2max` global variable in the `libgrasprci`
    !! module.
    function libgraspci_global_j2max() bind(c)
        integer(c_int) :: libgraspci_global_j2max
        libgraspci_global_j2max = j2max
    end function libgraspci_global_j2max

    !> Accessor function for the `NW` global variable in the `orb_C` module.
    function libgraspci_global_orb_nw() bind(c)
        use orb_C, only: NW
        integer(c_int) :: libgraspci_global_orb_nw
        libgraspci_global_orb_nw = NW
    end function libgraspci_global_orb_nw

    !> Accessor function for the `NCF` global variable in the `orb_C` module.
    function libgrasprci_global_orb_ncf() bind(c)
        use orb_C, only: NCF
        integer(c_int) :: libgrasprci_global_orb_ncf
        libgrasprci_global_orb_ncf = NCF
    end function libgrasprci_global_orb_ncf

    !> Accessor function for the `NVEC` global variable in the `prnt_C` module.
    function libgrasprci_global_prnt_nvec() bind(c)
        use prnt_C, only: NVEC
        integer(c_int) :: libgrasprci_global_prnt_nvec
        libgrasprci_global_prnt_nvec = NVEC
    end function libgrasprci_global_prnt_nvec

    !> Returns the `n` first values in `NP` and `NAK` arrays in `orb_C`.
    !!
    !! These arrays contain the contain the principal and \f$\kappa\f$ quantum
    !! numbers of the orbitals, respectively.
    !!
    !! @param n Length of the `np_c` and `nak_c` arrays.
    !! @param np_c Pointer to an array of at least length `n` which will be
    !!     filled with the first `n` principal quantum numbers.
    !! @param nak_c Pointer to an array of at least length `n` which will be
    !!     filled with the first `n` \f$\kappa\f$ quantum numbers.
    !!
    !! If the global `NW` is less than `n`, only the first `NW` values will be
    !! filled.
    subroutine libgraspci_global_orbitals(n, np_c, nak_c) bind(c)
        use orb_C, only: NW, NP, NAK
        integer(c_int), value :: n
        integer(c_int), intent(out) :: np_c(NW), nak_c(NW)
        n = min(n, NW)
        np_c(1:n) = NP(1:n)
        nak_c(1:n) = NAK(1:n)
    end

    !> Returns the C string of the last error
    !!
    !! Note: calling this function again will invalidate the previously returned
    !! pointer.
    function libgrasprci_error_string() bind(c)
        type(c_ptr) :: libgrasprci_error_string
        lasterror_cstr = trim(lasterror)//c_null_char
        libgrasprci_error_string = c_loc(lasterror_cstr)
    end

    !> Enable or disable debug output from `libgrasp-rci` routines.
    !!
    !! @param setting `true`/`false` enables/disables the debug output.
    subroutine libgrasprci_debug(setting) bind(c)
        logical(c_bool), value :: setting
        debug = setting
    end

    !===========================================================================
    ! The C API for working with the `matrix1pscalar` type
    !---------------------------------------------------------------------------

    !> Returns the `(i, j)` matrix element of a 1-particle scalar operator.
    !!
    !! @param matrix1pscalar_cptr The C pointer pointing to a `matrix1pscalar`
    !! object.
    !! @param i The row index of the matrix.
    !! @param j The column index of the matrix.
    !!
    !! _Note: no bounds checking is performed._
    function libgrasprci_matrix1pscalar(matrix1pscalar_cptr, i, j) bind(c)
        type(c_ptr), intent(in), value :: matrix1pscalar_cptr
        integer(c_int), intent(in), value :: i, j
        real(c_double) :: libgrasprci_matrix1pscalar
        type(matrix1pscalar), pointer :: matrix
        call c_f_pointer(matrix1pscalar_cptr, matrix)
        libgrasprci_matrix1pscalar = matrix%matrix(i, j)
    end function libgrasprci_matrix1pscalar

    !> Returns the number of orbitals in the basis that were used to calculate
    !! the matrix.
    !!
    !! @param matrix1pscalar_cptr The C pointer pointing to a `matrix1pscalar`
    !! object.
    !! @return The number of orbital (or, equivalently, the number of
    !! rows/columns of matrix).
    !!
    !! _Note: no bounds checking is performed._
    function libgrasprci_matrix1pscalar_n(matrix1pscalar_cptr) bind(c)
        type(c_ptr), intent(in), value :: matrix1pscalar_cptr
        integer(c_int) :: libgrasprci_matrix1pscalar_n
        type(matrix1pscalar), pointer :: matrix
        call c_f_pointer(matrix1pscalar_cptr, matrix)
        libgrasprci_matrix1pscalar_n = matrix%n
    end function libgrasprci_matrix1pscalar_n

    !> Deallocates a `matrix1pscalar` object referenced by `matrix1pscalar_cptr`.
    subroutine libgrasprci_matrix1pscalar_delete(matrix1pscalar_cptr) bind(c)
        type(c_ptr), intent(in), value :: matrix1pscalar_cptr
        type(matrix1pscalar), pointer :: matrix
        if(debug) then
            print '(a,z16.16)', "Deallocating matrix1pscalar: 0x", matrix1pscalar_cptr
        endif
        call c_f_pointer(matrix1pscalar_cptr, matrix)
        deallocate(matrix%matrix)
        deallocate(matrix)
    end subroutine libgrasprci_matrix1pscalar_delete

    !> Creates a copy of the `matrix1pscalar` referenced by `matrix1pscalar_cptr`
    !! C pointer.
    function libgrasprci_matrix1pscalar_copy(matrix1pscalar_cptr) bind(c)
        type(c_ptr), intent(in), value :: matrix1pscalar_cptr
        type(c_ptr) :: libgrasprci_matrix1pscalar_copy
        type(matrix1pscalar), pointer :: matrix, matrix_copy

        call c_f_pointer(matrix1pscalar_cptr, matrix)
        allocate(matrix_copy)
        allocate(matrix_copy%matrix(matrix%n,matrix%n))

        matrix_copy%n = matrix%n
        matrix_copy%matrix(:,:) = matrix%matrix(:,:)

        libgrasprci_matrix1pscalar_copy = c_loc(matrix_copy)
    end function libgrasprci_matrix1pscalar_copy

    !> Disables a contribution for the `i`th orbital by settings the
    !! corresponding row and column to zero.
    subroutine libgrasprci_matrix1pscalar_disable(matrix1pscalar_cptr, i) bind(c)
        use grasp_rciqed_kinds, only: dp

        type(c_ptr), intent(in), value :: matrix1pscalar_cptr
        integer(c_int), intent(in), value :: i

        integer :: k
        type(matrix1pscalar), pointer :: matrix

        call c_f_pointer(matrix1pscalar_cptr, matrix)
        do k = 1, matrix%n
            matrix%matrix(i, k) = 0.0_dp
            matrix%matrix(k, i) = 0.0_dp
        enddo
    end subroutine libgrasprci_matrix1pscalar_disable

    !===========================================================================
    ! Loading input files
    !---------------------------------------------------------------------------

    function libgrasprci_initalize(isodata, csl, orbitals, mixing) bind(c)
        use, intrinsic :: iso_c_binding
        use grasp_c_interface
        use grasp_rciqed_lib9290_init, only: lib9290_init_constants, lib9290_init_rkco_gg
        use grasp_rciqed_lib9290_files, only: load_isodata, load_csl, load_orbitals, load_mixing
        use grasp_rciqed, only: init_rkintc
        implicit none

        character(kind=c_char), intent(in) :: isodata(1), csl(1), orbitals(1), mixing(1)
        integer(c_int) :: libgrasprci_initalize

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

        libgrasprci_initalize = 0
    end function libgrasprci_initalize

    !===========================================================================
    ! C API for accessing the orbitals, ASFs coefficients etc.
    !---------------------------------------------------------------------------

    !> Returns the coefficient of the `i`th CSF in the `k`th ASF.
    !!
    !! @param k Index of the ASF (`1 <= k <= NVEC`)
    !! @param i Index of the CSF (many-body basis state; `1 <= i <= NCF`)
    !! @return The corresponding expansion coefficient of the specified ASF.
    !!
    !! Internally, it accesses the appropriate element on the `eigv_C::EVEC`
    !! array.
    function libgrasprci_asfcoefficient(k, i) bind(c)
        use eigv_C, only: EVEC
        use orb_C, only: NCF
        integer(c_int), intent(in), value :: k, i
        real(c_double) :: libgrasprci_asfcoefficient
        libgrasprci_asfcoefficient = EVEC(i + (k - 1) * NCF)
    end function libgrasprci_asfcoefficient

    !===========================================================================
    ! C API for calculating various physical operators
    !---------------------------------------------------------------------------

    !> Allocates and stores a C pointer to a `matrix1pscalar` type in
    !! `matrix_cptr`, filled with the 1-particle matrix elements of Dirac
    !! operator + the central potential.
    !!
    !! @param matrix_cptr A reference to a pointer that will be set to the
    !!     address of the `matrix1pscalar` object.
    !! @return Non-zero if there was an error.
    function libgrasprci_diracpot_matrix1pscalar(matrix_cptr) bind(c)
        use grasp_rciqed_kinds, only: dp
        use orb_C, only: NW, NAK
        type(c_ptr), intent(out) :: matrix_cptr
        integer(c_int) :: libgrasprci_diracpot_matrix1pscalar
        type(matrix1pscalar), pointer :: matrix
        integer :: i, j

        allocate(matrix)
        matrix%n = NW
        allocate(matrix%matrix(NW,NW))
        do i = 1, NW
            do j = 1, NW ! TODO: only calculate one half of the matrix
                if(NAK(i) == NAK(j)) then
                    call iabint_safe(i, j, matrix%matrix(i, j))
                else
                    matrix%matrix(i, j) = 0.0_dp
                endif
            enddo
        enddo
        matrix_cptr = c_loc(matrix)
        libgrasprci_diracpot_matrix1pscalar = 0 ! no error
    end function libgrasprci_diracpot_matrix1pscalar

    !> Calls the `IABINT` routine, but ensures that the input `k` and `l`
    !! arguments do not get swapped.
    subroutine iabint_safe(k, l, result)
        use grasp_rciqed_kinds, only: real64
        use iabint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call IABINT(k, l, result)
    end

    !> Allocates and stores a C pointer to a `matrix1pscalar` type in
    !! `matrix_cptr`, filled with the 1-particle matrix elements of the QED
    !! self-energy operators corresponding to `setype`.
    !!
    !! @param setype Integer specifying the self-energy type.
    !! @param matrix_cptr A reference to a pointer that will be set to the
    !!     address of the `matrix1pscalar` object.
    !! @return Non-zero if there was an error.
    function libgrasprci_qedse_matrix1pscalar(setype, matrix_cptr) bind(c)
        use orb_C, only: NW
        use grasp_rciqed_qed, only: qedse

        integer(c_int), intent(in), value :: setype
        type(c_ptr), intent(out) :: matrix_cptr
        integer(c_int) :: libgrasprci_qedse_matrix1pscalar

        type(matrix1pscalar), pointer :: sematrix

        if(setype < 1 .or. setype > 4) then
            lasterror = "grasp_ci_qedse_matrix: invalid setype specified"
            libgrasprci_qedse_matrix1pscalar = 1
            return
        endif

        allocate(sematrix)
        sematrix%n = NW
        allocate(sematrix%matrix(NW,NW))
        call qedse(setype, sematrix%matrix)
        matrix_cptr = c_loc(sematrix)

        libgrasprci_qedse_matrix1pscalar = 0 ! no error
    end function libgrasprci_qedse_matrix1pscalar

    !===========================================================================
    ! C API for working with the onescalar_cache type
    !---------------------------------------------------------------------------

    !> Run the `ONESCALAR` routine for the CSF pair `(ir, ic)` and return the
    !! output values as a `onescalar_cache` object.
    !!
    !! @param ir Row index of a CSF.
    !! @param ic Column index of a CSF.
    !! @return A C pointer to an `onescalar_cache` object.
    function libgrasprci_onescalar(ir, ic) bind(c)
        use orb_C, only: NW
        use onescalar_I
        integer(c_int), intent(in), value :: ir, ic
        type(c_ptr) :: libgrasprci_onescalar
        type(onescalar_cache), pointer :: cache

        allocate(cache)
        allocate(cache%tshell(NW))
        cache%ic = ic
        cache%ir = ir
        call ONESCALAR(ic, ir, cache%k, cache%l, cache%tshell)
        libgrasprci_onescalar = c_loc(cache)
    end

    !> Deallocates the `onescalar_cache` object referenced by `onescalar_cptr`.
    subroutine libgrasprci_onescalar_delete(onescalar_cptr) bind(c)
        type(c_ptr), intent(in), value :: onescalar_cptr
        type(onescalar_cache), pointer :: onescalar
        if(debug) then
            print '(a,z16.16)', "Deallocating onescalar_cache: 0x", onescalar_cptr
        endif
        call c_f_pointer(onescalar_cptr, onescalar)
        deallocate(onescalar%tshell)
        deallocate(onescalar)
    end subroutine libgrasprci_onescalar_delete

    !===========================================================================
    ! C API for calculating matrix elements
    !---------------------------------------------------------------------------

    function libgrasprci_matrixelement_1p(cache_cptr, matrix_cptr) bind(c)
        use grasp_rciqed_cimatrixelements, only: qed_se
        type(c_ptr), value :: cache_cptr, matrix_cptr
        real(c_double) :: libgrasprci_matrixelement_1p
        type(onescalar_cache), pointer :: cache
        type(matrix1pscalar), pointer :: matrix

        call c_f_pointer(cache_cptr, cache)
        call c_f_pointer(matrix_cptr, matrix)
        libgrasprci_matrixelement_1p = qed_se(matrix%matrix, cache%ic, cache%ir, cache%k, cache%l, cache%tshell)
    end function libgrasprci_matrixelement_1p

    function libgrasprci_matrixelement_diracpot(cache_cptr) bind(c)
        use grasp_rciqed_cimatrixelements, only: dirac_potential, coulomb
        type(c_ptr), value :: cache_cptr
        real(c_double) :: libgrasprci_matrixelement_diracpot
        type(onescalar_cache), pointer :: cache
        call c_f_pointer(cache_cptr, cache)
        libgrasprci_matrixelement_diracpot = dirac_potential(cache%ic, cache%ir, cache%k, cache%l, cache%tshell)
    end function libgrasprci_matrixelement_diracpot

    function libgrasprci_matrixelement_coulomb(ir, ic) bind(c)
        use grasp_rciqed_cimatrixelements, only: coulomb
        integer(c_int), value :: ir, ic
        real(c_double) :: libgrasprci_matrixelement_coulomb
        libgrasprci_matrixelement_coulomb = coulomb(ic, ir)
    end function libgrasprci_matrixelement_coulomb

end module libgrasprci
