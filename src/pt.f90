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
    use grasp_rciqed_qed, only: init_vacuum_polarization, qedse, &
        nsetypes, setypes_long, setypes_short
    use grasp_rciqed_rcisettings, only: rcisettings, read_settings_toml
    use grasp_rciqed_cimatrixelements
    use decide_C, only: LSE
    use orb_C, only: NW
    use prnt_C, only: NVEC
    use qedcut_C, only: NQEDCUT, NQEDMAX
    implicit none

    integer :: i
    type matrixelement
        real(real64) :: diracpot = 0.0_dp, coulomb = 0.0_dp, breit = 0.0_dp, &
            nms = 0.0_dp, sms = 0.0_dp, vp = 0.0_dp, se = 0.0_dp, &
            se_array(nsetypes) = (/(0.0_dp, i = 1, nsetypes)/)
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

    character(:), allocatable :: file_csls, file_wfns, file_mixing, state
    type(rcisettings) :: settings

    integer :: status, k, j, j2max
    integer :: fh

    type(matrixelement) :: hij, asfv
    type(matrixelement), dimension(:), allocatable :: asfvalues
    real(real64), allocatable :: sematrix(:,:,:)
    real(real64) :: sum

    type(integer), dimension(:), allocatable :: hcs_rci, hcs_pt
    type(logical), dimension(hcs_len) :: hcs_cat

    character(:), allocatable :: cli_argument, cli_cutoff_value_string
    logical :: cli_status
    logical :: cli_state_found = .false.
    character(*), parameter :: cli_cutoff_name = '--qed-se-hydrogenic-cutoff='
    integer :: cli_cutoff_value = -1

    character(*), parameter :: isodata = 'isodata' ! name of the isodata file

    ! Parse the command line arguments
    do i = 1, command_argument_count()
        cli_status = get_command_argument_allocating(i, cli_argument)
        if(.not.cli_status) then
            error stop 'Error while parsing command line arguments'
        endif
        if(len(cli_argument) > 0 .and. cli_argument(1:1) == "-") then
            if(cli_argument(1:len(cli_cutoff_name)) == cli_cutoff_name) then
                cli_cutoff_value_string = cli_argument(len(cli_cutoff_name)+1:len(cli_argument))
                read(cli_cutoff_value_string, *, iostat=status) cli_cutoff_value
                if(status /= 0) then
                    print '(a," ",a)', "ERROR: Parsing an option failed:", cli_argument
                    print '(a," ",a)', "Unable to convert to integer:", cli_cutoff_value_string
                    error stop 'Error while parsing command line arguments'
                endif
                if(cli_cutoff_value < 0) then
                    print '(a," ",a)', "ERROR: Parsing an option failed:", cli_argument
                    print '(a," ",i0)', "Cutoff must be >= 0. Passed:", cli_cutoff_value
                    error stop 'Error while parsing command line arguments'
                endif
            else
                print '(a," ",a)', "ERROR: Invalid option passed:", cli_argument
                error stop 'Error while parsing command line arguments'
            endif
        elseif(.not.cli_state_found) then
            state = cli_argument
            cli_state_found = .true.
        else
            print '(a)', "ERROR: Multiple positional arguments passed."
            error stop 'Error while parsing command line arguments'
        endif
    enddo
    if(.not.cli_state_found) then
        print '(a)', "ERROR: STATE positional argument not passed."
        error stop 'Error while parsing command line arguments'
    endif

    ! The input files are assumed to be <state>.c, <state>.w and <state>.cm
    file_csls = trim(state)//'.c'
    file_wfns = trim(state)//'.w'
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

    ! Optionally, override the QED SE cutoff value with command line arguments
    if(cli_cutoff_value >= 0 .and. settings%qed_se_hydrogenic_cutoff >= 0) then
        print '(a)', "WARNING: Hydrogenic SE cutoff overriden with command line argument"
    elseif(cli_cutoff_value >= 0 .and. settings%qed_se_hydrogenic_cutoff < 0) then
        settings%qed_se_hydrogenic_cutoff = cli_cutoff_value
    endif

    print "(a)", "Hamiltonian parts enabled in RCI calculation:"
    print '(" ",a,"  ",a)', check(.true.), "Dirac + nuclear potential"
    print '(" ",a,"  ",a)', check(settings%breit_enabled), "Breit"
    print '(" ",a,"  ",a)', check(settings%nms_enabled), "Normal mass shift"
    print '(" ",a,"  ",a)', check(settings%sms_enabled), "Special mass shift"
    print '(" ",a,"  ",a)', check(settings%qed_vp_enabled), "QED vacuum polarization"
    if(settings%qed_se /= 0) then
        print '(" ",a,"  ",a,": ",a," (#",i0,": ",a,")")', &
            check(.true.), "QED self-energy", trim(setypes_short(settings%qed_se)), &
            settings%qed_se, trim(setypes_long(settings%qed_se))
    else
        print '(" ",a,"  ",a)', check(.false.), "QED self-energy"
    endif
    if(settings%qed_se_hydrogenic_cutoff >= 0) then
        print '("    ",a,":  ",i0)', "Hydrogenic SE n cutoff", settings%qed_se_hydrogenic_cutoff
    endif

    hcs_cat(1) = .true.
    hcs_cat(2) = .true.
    hcs_cat(3) = settings%breit_enabled
    hcs_cat(4) = settings%nms_enabled
    hcs_cat(5) = settings%sms_enabled
    hcs_cat(6) = settings%qed_vp_enabled
    hcs_cat(7) = .not.(settings%qed_se == 0)

    ! Also set the hydrogenic settings in qedcut_C appropriately:
    ! LSE needs to be set to .true. since QED_SLFEN checks for it when applying
    ! NQEDMAX for some reason.
    LSE = .true.
    if(settings%qed_se_hydrogenic_cutoff >= 0) then
        NQEDCUT = 1
        NQEDMAX = settings%qed_se_hydrogenic_cutoff
    else
        NQEDCUT = 0
    endif

    ! Calculate the self-energy matrix (or matrices, if SE in in PT mode).
    ! If settings%qed_se == 0, then we need to calculate _all_ SE operators.
    if(settings%qed_se == 0) then
        allocate(sematrix(nsetypes,NW,NW))
        do i = 1, nsetypes
            call qedse(i, sematrix(i,:,:))
        enddo
    else
        allocate(sematrix(1,NW,NW))
        call qedse(settings%qed_se, sematrix(1,:,:))
    endif

    ! Evaluate the ASF expectation values and write them to a CSV file.
    print '(a)', ">>> Evaluating ASF expectation values"
    allocate(asfvalues(NVEC))

    do i = 1, ncsfs_global()
        do j = 1, i
            ! Evaluate the (i, j) CI matrix elements for all the components of
            ! the Hamiltonian
            hij%diracpot = dirac_potential(i, j)
            hij%coulomb = coulomb(i, j)
            hij%breit = breit(i, j)
            hij%nms = nms(i, j)
            hij%sms = sms(i, j)
            hij%vp = qed_vp(i, j)

            if(settings%qed_se == 0) then
                do k = 1, nsetypes
                    hij%se_array(k) = qed_se(sematrix(k,:,:), i, j)
                enddo
            else
                hij%se = qed_se(sematrix(1,:,:), i, j)
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
        print '(a1,a16,a18)',    " ", "Type", "<H_i> (Ha)"
        sum = 0.0_dp
        do i = 1, hcs_len
            if(.not.hcs_cat(i)) cycle
            if(i == 7) then
                ! Self-energy gets special treatment
                print '("QED SE ",a10,1p,e18.8)', &
                    trim(setypes_short(settings%qed_se)), getindex(asfv, i)
            else
                print '(a1,a16,1p,e18.8)', " ", trim(hcs_pretty(i)), getindex(asfv, i)
            endif
            sum = sum + getindex(asfv, i)
        enddo
        print '(a1,a16,1p,e18.8)', ">", "DC", asfv%diracpot + asfv%coulomb
        print '(a1,a16,1p,e18.8)', ">", "Non.pet", sum
        do i = 1, hcs_len
            if(hcs_cat(i)) cycle
            if(i == 7) then
                ! Self-energy gets special treatment, since in PT mode, we have
                ! several parallel methods.
                do j = 1, nsetypes
                    print '("QED SE ",a10,1p,e18.8)', trim(setypes_short(j)), asfv%se_array(j)
                enddo
            else
                print '(a1,a16,1p,e18.8)', " ", trim(hcs_pretty(i)), getindex(asfv, i)
            endif
        enddo
    enddo

    close(fh)

contains

    function getindex(me, idx)
        use grasp_rciqed_kinds, only: real64
        type(matrixelement), intent(in) :: me
        integer, intent(in) :: idx
        real(real64) :: getindex

        if(idx < 1 .or. idx > 8) then
            error stop "Bad index in getindex(::matrixelement)"
        endif

        if(idx == 1) getindex = me%diracpot
        if(idx == 2) getindex = me%coulomb
        if(idx == 3) getindex = me%breit
        if(idx == 4) getindex = me%nms
        if(idx == 5) getindex = me%sms
        if(idx == 6) getindex = me%vp
        if(idx == 7) getindex = me%se
    end function getindex

    subroutine write_csv_header(unit)
        integer, intent(in) :: unit

        write(unit, '(a9)', advance='no') "asf_index"
        do i = 1, hcs_len
            if(.not.hcs_cat(i)) cycle
            if(i == 7) then
                ! Self-energy gets special treatment
                write(unit, '(",",a24)', advance='no') trim(hcs_technical(i))//":"//trim(setypes_short(settings%qed_se))
            else
                write(unit, '(",",a24)', advance='no') trim(hcs_technical(i))
            endif
        enddo
        write(unit, '(",",a24)', advance='no') "sum:non_pt"
        do i = 1, hcs_len
            if(hcs_cat(i)) cycle
            if(i == 7) then
                ! Self-energy gets special treatment, since in PT mode, we have
                ! several parallel methods.
                do j = 1, nsetypes
                    write(unit, '(",",a24)', advance='no') "pt:"//trim(hcs_technical(i))//":"//trim(setypes_short(j))
                enddo
            else
                write(unit, '(",",a24)', advance='no') "pt:"//trim(hcs_technical(i))
            endif
        enddo
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
            if(i == 7) then
                ! Self-energy gets special treatment, since in PT mode, we have
                ! several parallel methods.
                do j = 1, nsetypes
                    write(unit, '(1(",",e24.16))', advance='no') me%se_array(j)
                enddo
            else
                write(unit, '(1(",",e24.16))', advance='no') getindex(me, i)
            endif
        enddo
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
        type(matrixelement), intent(in) :: hij
        type(matrixelement), intent(inout) :: asfvalue

        real(real64) :: cicj

        cicj = asf_coefficient(k, i) * asf_coefficient(k, j)
        asfvalue%diracpot = asfvalue%diracpot + hij%diracpot * cicj
        asfvalue%coulomb = asfvalue%coulomb + hij%coulomb * cicj
        asfvalue%breit = asfvalue%breit + hij%breit * cicj
        asfvalue%vp  = asfvalue%vp  + hij%vp * cicj
        asfvalue%nms = asfvalue%nms + hij%nms * cicj
        asfvalue%sms = asfvalue%sms + hij%sms * cicj
        asfvalue%se = asfvalue%se + hij%se * cicj
        asfvalue%se_array(:) = asfvalue%se_array(:) + hij%se_array(:) * cicj
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

        ! If the values are exactly the same, we also say that they are within tolerance.
        if (abs(a-b) == 0) then
            within_tolerance = .true.
            return
        endif

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

    !> Retrieves the `n`-th command line argument.
    !!
    !! If it was able to fetch the argument value, returns `.true.` and sets
    !! `value` to the value ([re]allocating if necessary). Returns `.false.` if
    !! there was an error retrieving the argument.
    !!
    !! Under the hood it calls `get_command_argument`, but it properly allocates
    !! or re-allocates the `value` to match the actual length of the argument.
    !!
    !! TODO: For this to be a proper library function, it should be implemented
    !! as an interface with additional methods to handle pointers and fixed-length
    !! strings.
    function get_command_argument_allocating(n, value)
        logical :: get_command_argument_allocating
        integer, intent(in) :: n
        character(:), allocatable, intent(inout) :: value
        integer :: length, status
        character(1) :: test

        ! First, make an inquiry call to get_environment_variable to determine
        ! whether the variable exists and, if so, its length.
        call get_command_argument(n, test, length, status)
        if(status > 0) then
            ! From the GFortran manual:
            !
            ! > If the argument retrieval fails, STATUS is a positive number;
            !
            ! So a positive status will make this function fail (i.e. return
            ! .false.)
            get_command_argument_allocating = .false.
            return
        endif
        ! Will allocate or re-allocate value, unless it already has the correct length.
        if(allocated(value) .and. len(value) /= length) deallocate(value)
        if(.not.allocated(value)) allocate(character(length) :: value)
        ! Only call get_command_argument again if the argument is non-empty
        if(length > 0) then
            call get_command_argument(n, value, length, status)
        endif
        get_command_argument_allocating = .true.
    end

end program rci_qed_pt
