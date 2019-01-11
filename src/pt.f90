!>     ./rci-qed.pt <state>
!!
!! Performs analysis and perturbative calculation on the RCI state.
!!
!! Takes a single command line argument: `state`.
!!
!! Depends on files: `isodata`, `<state>.c`, `<state>.w` and `<state>.cm`
!!
program rci_qed_pt
    use grasp_rciqed_kinds, only: real64, dp
    use grasp_lib9290, only: init_isocw_full
    use grasp_rciqed_lib9290_files, only: load_mixing
    use grasp_rciqed_lib9290_csls, only: ncsfs_global
    use grasp_rciqed, only: init_rkintc
    use grasp_rciqed_breit, only: init_breit
    use grasp_rciqed_mass_shifts, only: init_mass_shifts
    use grasp_rciqed_qed, only: init_vacuum_polarization
    use grasp_rciqed_rcisettings, only: rcisettings, read_settings_toml
    use grasp_cimatrixelements
    use prnt_C, only: NVEC
    implicit none

    type matrixelement
        real(real64) :: diracpot, coulomb, breit, nms, sms, vp, se_mohr
    end type matrixelement

    integer, parameter :: hcs_len = 7
    character(*), dimension(hcs_len), parameter :: hcs_pretty = (/ &
        "KE + Nuc.pot", "Coulomb     ", "Breit       ", "NMS         ", "SMS         ", &
        "QED VP      ", "QED SE      " &
    /)
    character(*), dimension(hcs_len), parameter :: hcs_technical = (/ &
        "kinetic_nucl", "coulomb     ", "breit       ", "nms         ", "sms         ", &
        "qed_vp      ", "qed_se      " &
    /)

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls, file_wfns, file_mixing
    type(rcisettings) :: settings

    integer :: status, k, i, j, j2max
    integer :: fh

    type(matrixelement) :: hij, asfv
    type(matrixelement), dimension(:), allocatable :: asfvalues
    real(real64) :: sum

    type(integer), dimension(:), allocatable :: hcs_rci, hcs_pt
    type(logical), dimension(hcs_len) :: hcs_cat

    character(*), parameter :: isodata = 'isodata' ! name of the isodata file

    ! Fetch the first command line argument, which should be the name of the state.
    if (command_argument_count() /= 1) then
        error stop 'Missing command line argument: ./rci_qedreport <STATE>'
    endif
    call get_command_argument(1, state) ! TODO: make state allocatable
    state_len = len_trim(state)

    ! The input files are assumed to be <state>.c, <state>.w and <state>.cm
    allocate(character(state_len+2) :: file_csls)
    file_csls = trim(state)//'.c'

    allocate(character(state_len+2) :: file_wfns)
    file_wfns = trim(state)//'.w'

    ! TODO: should be +3? Or rather, allocate is not actually necessary?
    allocate(character(state_len+2) :: file_mixing)
    file_mixing = trim(state)//'.cm'

    ! Make sure that all the files exist, or ERROR STOP otherwise
    call check_file(isodata)
    call check_file(file_csls)
    call check_file(file_wfns)
    call check_file(file_mixing)
    call check_file(trim(state)//".settings.toml")

    ! Load the isotope information, CSL, radial orbitals and the CI mixing coefficents.
    print '(a)', ">>> CALLING: grasp_lib9290::init_isocw_full"
    call init_isocw_full(isodata, file_csls, file_wfns)
    print '(a)', ">>> CALLING: lib92_init_mixing"
    call load_mixing(file_mixing)
    print '(a)', ">>> INITIALIZING rci commons"
    call init_rkintc(j2max)
    call init_breit(j2max)
    call init_vacuum_polarization
    call init_mass_shifts

    ! Load the info in $(state).settings.toml file
    print '(">>> LOADING: ",a,".settings.toml")', trim(state)
    if(.not.read_settings_toml(trim(state), settings)) then
        print '("ERROR: failed to read ",a,".setting.toml")', trim(state)
        error stop
    endif
    call verify_rcisettings(settings)

    print "(a)", "Hamiltonian parts enabled in RCI calculation:"
    print '(" ",a,"  ",a)', check(.true.), "Dirac + nuclear potential"
    print '(" ",a,"  ",a)', check(settings%breit_enabled), "Breit"
    print '(" ",a,"  ",a)', check(settings%nms_enabled), "Normal mass shift"
    print '(" ",a,"  ",a)', check(settings%sms_enabled), "Special mass shift"
    print '(" ",a,"  ",a)', check(settings%qed_vp_enabled), "QED vacuum polarization"
    print '(" ",a,"  ",a)', check(settings%qed_se_enabled), "QED self-energy"

    hcs_cat(1) = .true.
    hcs_cat(2) = .true.
    hcs_cat(3) = settings%breit_enabled
    hcs_cat(4) = settings%nms_enabled
    hcs_cat(5) = settings%sms_enabled
    hcs_cat(6) = settings%qed_vp_enabled
    hcs_cat(7) = settings%qed_se_enabled

    ! Evaluate the ASF expectation values and write them to a CSV file.
    print '(a)', ">>> Evaluating ASF expectation values"
    allocate(asfvalues(NVEC))
    do k = 1, NVEC
        call matrixelement_zero(asfvalues(k))
    enddo

    do i = 1, ncsfs_global()
        do j = 1, i
            ! Evaluate the (i, j) CI matrix elements for all the components of
            ! the Hamiltonian
            hij%diracpot = dirac_potential(i, j)
            hij%coulomb = coulomb(i, j)
            hij%breit = breit(i, j)
            hij%vp = qed_vp(i, j)
            hij%nms = nms(i, j)
            hij%sms = sms(i, j)
            if (i == j) then
                hij%se_mohr = qed_se_mohr(i)
            else
                hij%se_mohr = 0.0_dp
            endif

            ! Add the values to the ASF values
            do k = 1, NVEC
                call add_asfvalue(i, j, hij, k, asfvalues(k))
                if(i /= j) call add_asfvalue(j, i, hij, k, asfvalues(k))
            enddo
        enddo
    enddo
    close(fh)

    print '(a)', ">>> Outputting data:"
    open(newunit=fh, file=trim(state)//".asfenergies.csv", action='write')
    call write_csv_header(fh)

    do k = 1, NVEC
        asfv = asfvalues(k)

        ! Write to the CSV file
        call write_csv_asf_line(fh, asfv)

        ! Also print them out
        print '(80("="))'
        print '(">>> ASF #",i0)', k
        print '(80("-"))'
        print '(a1,a12,a18)',    " ", "Type", "<H_i> (Ha)"
        sum = 0.0_dp
        do i = 1, hcs_len
            if(.not.hcs_cat(i)) cycle
            print '(a1,a12,1p,e18.8)', " ", trim(hcs_pretty(i)), getindex(asfv, i)
            sum = sum + getindex(asfv, i)
        enddo
        print '(a1,a12,1p,e18.8)', ">", "DC", asfv%diracpot + asfv%coulomb
        print '(a1,a12,1p,e18.8)', ">", "Non.pet", sum
        do i = 1, hcs_len
            if(hcs_cat(i)) cycle
            print '(a1,a12,1p,e18.8)', " ", trim(hcs_pretty(i)), getindex(asfv, i)
            sum = sum + getindex(asfv, i)
        enddo
        print '(a1,a12,1p,e18.8)', ">", "Full", &
            asfv%diracpot + asfv%coulomb + asfv%breit + asfv%vp + asfv%nms + asfv%sms + asfv%se_mohr
    enddo

    close(fh)

contains

    function getindex(me, idx)
        use grasp_rciqed_kinds, only: real64
        type(matrixelement), intent(in) :: me
        integer, intent(in) :: idx
        real(real64) :: getindex

        if(idx < 1 .or. idx > 7) then
            error stop "Bad index in getindex(::matrixelement)"
        endif

        if(idx == 1) getindex = me%diracpot
        if(idx == 2) getindex = me%coulomb
        if(idx == 3) getindex = me%breit
        if(idx == 4) getindex = me%nms
        if(idx == 5) getindex = me%sms
        if(idx == 6) getindex = me%vp
        if(idx == 7) getindex = me%se_mohr
    end function getindex

    subroutine write_csv_header(unit)
        integer, intent(in) :: unit

        write(unit, '(a9)', advance='no') "asf_index"
        do i = 1, hcs_len
            if(.not.hcs_cat(i)) cycle
            write(unit, '(",",a24)', advance='no') trim(hcs_technical(i))
        enddo
        write(unit, '(",",a24)', advance='no') "sum:non_pt"
        do i = 1, hcs_len
            if(hcs_cat(i)) cycle
            write(unit, '(",",a24)', advance='no') "pt:"//trim(hcs_technical(i))
        enddo
        write(unit, '(",",a24)', advance='no') "sum:with_pt"
        write(unit, '()')
    end subroutine write_csv_header

    subroutine write_csv_asf_line(unit, me)
        integer, intent(in) :: unit
        type(matrixelement), intent(in) :: me

        real(real64) :: sum = 0.0_dp

        write(unit, '(i9)', advance='no') k
        do i = 1, hcs_len
            if(.not.hcs_cat(i)) cycle
            write(unit, '(1(",",e24.16))', advance='no') getindex(me, i)
            sum = sum + getindex(me, i)
        enddo
        write(unit, '(1(",",e24.16))', advance='no') sum
        do i = 1, hcs_len
            if(hcs_cat(i)) cycle
            write(unit, '(1(",",e24.16))', advance='no') getindex(me, i)
            sum = sum + getindex(me, i)
        enddo
        write(unit, '(1(",",e24.16))', advance='no') sum
        write(unit, '()') ! write the newline
    end subroutine write_csv_asf_line

    subroutine check_file(filename)
        use grasp_rciqed_system, only: file_exists
        character(*), intent(in) :: filename

        if(.not.file_exists(filename)) then
            print '("Missing file ",a)', filename
            error stop
        endif
    end subroutine check_file

    !> Zeroes all the fields of a `matrixelement`.
    subroutine matrixelement_zero(h)
        use grasp_rciqed_kinds, only: dp

        type(matrixelement), intent(out) :: h

        h%diracpot = 0.0_dp
        h%coulomb = 0.0_dp
        h%breit = 0.0_dp
        h%vp = 0.0_dp
        h%nms = 0.0_dp
        h%sms = 0.0_dp
        h%se_mohr = 0.0_dp
    end subroutine matrixelement_zero

    !> Adds the contribution of the `(i, j)` matrix element to the `asfvalue`.
    !!
    !! `hij` must contain the `(i, j)` Hamiltonian matrix element and asfvalue
    !! must be the `k`-th ASF value.
    subroutine add_asfvalue(i, j, hij, k, asfvalue)
        use grasp_rciqed_kinds, only: real64, dp
        use orb_C
        use prnt_C
        use syma_C
        use eigv_C
        use eigvec1_C
        implicit none

        integer, intent(in) :: i, j, k
        type(matrixelement), intent(out) :: hij
        type(matrixelement), intent(inout) :: asfvalue

        real(real64) :: cicj

        cicj = asf_coefficient(k, i) * asf_coefficient(k, j)
        asfvalue%diracpot = asfvalue%diracpot + hij%diracpot * cicj
        asfvalue%coulomb = asfvalue%coulomb + hij%coulomb * cicj
        asfvalue%breit = asfvalue%breit + hij%breit * cicj
        asfvalue%vp  = asfvalue%vp  + hij%vp * cicj
        asfvalue%nms = asfvalue%nms + hij%nms * cicj
        asfvalue%sms = asfvalue%sms + hij%sms * cicj
        asfvalue%se_mohr = asfvalue%se_mohr + hij%se_mohr * cicj
    end subroutine add_asfvalue

    !> Returns the i-th ASF coefficient of the k-th state
    function asf_coefficient(k, i)
        use eigv_C, only: EVEC
        use orb_C, only: NCF
        integer, intent(in) :: k, i
        real(real64) :: asf_coefficient
        asf_coefficient = EVEC(i + (k - 1) * NCF)
    end function asf_coefficient

    !> Verifies if the nuclear and grid values in the `.settings.toml` file match
    !! the ones set by loading the `isodata` files.
    !!
    !! Does not abort if there are discrepancies, just prints warnings.
    subroutine verify_rcisettings(settings)
        use grasp_rciqed_kinds, only: real64, dp
        use grasp_rciqed_rcisettings, only: rcisettings
        use decide_C, only: LTRANS, LNMS, LSMS, LVP, LSE
        use def_C, only: Z, EMN, AUMAMU, FMTOAU
        use grid_C, only: RNT, H, N
        use npar_C, only: NPARM, PARM

        real(real64), parameter :: rtol = 1e-15_dp

        type(rcisettings), intent(in) :: settings

        if(.not.within_tolerance(settings%Z, Z, rtol)) then
            print '(a)', "WARNING: discrepancy between isodata and .settings.toml"
            print '(a)', "  values for Z do not match"
            print '(a,e22.16)', "  Z in isodata (global): ", Z
            print '(a,e22.16)', "  Z in .settings.toml  : ", settings%Z
        endif

        if(.not.within_tolerance(settings%atomic_mass_amu, EMN * AUMAMU, rtol)) then
            print '(a)', "WARNING: discrepancy between isodata and .settings.toml"
            print '(a)', "  values for atomic_mass_amu <-> EMN do not match"
            print '(a,e22.16)', "  EMN * AUMAMU (global)            : ", EMN * AUMAMU
            print '(a,e22.16)', "  atomic_mass_amu in .settings.toml: ", settings%atomic_mass_amu
        endif

        if(.not.within_tolerance(settings%atomic_mass_amu, EMN * AUMAMU, rtol)) then
            print '(a)', "WARNING: discrepancy between isodata and .settings.toml"
            print '(a)', "  values for atomic_mass_amu <-> EMN do not match"
            print '(a,e22.16)', "  EMN * AUMAMU (global)            : ", EMN * AUMAMU
            print '(a,e22.16)', "  atomic_mass_amu in .settings.toml: ", settings%atomic_mass_amu
        endif

    end subroutine verify_rcisettings

    !> Checks if the difference of `a` and `b` are within the tolerance relative
    !! to \f$\max(|a|,|b|)\f$.
    !!
    !! @param a,b Values to be checked.
    !! @param relative_tolerance Relative tolerance \f$\sigma\f$.
    !! @returns Whether \f$|a-b| / \max(|a|,|b|) < \sigma\f$.
    function within_tolerance(a, b, relative_tolerance)
        use grasp_rciqed_kinds, only: real64

        real(real64), intent(in) :: a, b, relative_tolerance
        logical :: within_tolerance
        real(real64) :: relative_difference

        relative_difference = abs(a-b) / max(abs(a), abs(b))
        if (relative_difference < relative_tolerance) then
            within_tolerance = .true.
        else
            within_tolerance = .false.
        endif
    end function within_tolerance

    !> Calculate the relative difference of `a` and `b`.
    !!
    !! The relative difference is defined as:
    !!
    !! \f[
    !!   \frac{|a-b|}{\max(|a|, |b|)}
    !! \f]
    !!
    !! @param a,b Input values.
    !! @returns The relative difference of `a` and `b`.
    function reldiff(a, b)
        use grasp_rciqed_kinds, only: real64
        real(real64), intent(in) :: a, b
        real(real64) :: reldiff
        reldiff = abs(a-b) / max(abs(a), abs(b))
    end function reldiff

    pure function check(yesno)
        logical, intent(in) :: yesno
        character(:), allocatable :: check
        if(yesno) then
            check = "✔"
        else
            check = "✘"
        endif
    end function check

end program rci_qed_pt
