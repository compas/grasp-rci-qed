!>     ./rci.matrixelements <state>
!!
!! Performs analysis and perturbative calculation on the RCI state.
!!
!! Takes a single command line argument: `state`.
!!
!! Depends on files: `isodata`, `<state>.c`, `<state>.w` and `<state>.cm`
!!
program matrixelements
    use grasp_kinds, only: real64, dp
    use grasp_lib9290, only: init_isocw_full
    use grasp_lib9290_files, only: load_mixing
    use grasp_lib9290_csls, only: ncsfs_global
    use grasp_rciqed, only: init_rkintc
    use grasp_rciqed_breit, only: init_breit
    use grasp_rciqed_mass_shifts, only: init_mass_shifts
    use grasp_rciqed_qed, only: init_vacuum_polarization
    use grasp_rciqed_rcisettings, only: rcisettings, read_settings_toml
    use grasp_cimatrixelements
    use prnt_C, only: NVEC
    implicit none

    type matrixelement
        real(real64) :: diracpot, coulomb, breit, vp, nms, sms, se_mohr
    end type matrixelement

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls, file_wfns, file_mixing
    type(rcisettings) :: settings

    integer :: status, k, i, j, j2max
    integer :: fh

    type(matrixelement) :: hij
    type(matrixelement), dimension(:,:), allocatable :: hamiltonian
    real(real64) :: result
    real(real64), dimension(20) :: tshell

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

    ! Allocate the dense CI matrix
    allocate(hamiltonian(ncsfs_global(), ncsfs_global()))

    print '(a)', ">>> Calculating matrixelements"
    do i = 1, ncsfs_global()
        do j = 1, i
            hamiltonian(i,j)%diracpot = dirac_potential(i, j)
            hamiltonian(i,j)%coulomb = coulomb(i, j)
            hamiltonian(i,j)%breit = breit(i, j)
            hamiltonian(i,j)%vp = qed_vp(i, j)
            hamiltonian(i,j)%nms = nms(i, j)
            hamiltonian(i,j)%sms = sms(i, j)
            if (i == j) then
                hamiltonian(i,j)%se_mohr = qed_se_mohr(i)
            else
                hamiltonian(i,j)%se_mohr = 0.0_dp
                hamiltonian(j,i) = hamiltonian(i,j)
            endif
        enddo
    enddo
    close(fh)

    ! Evaluate the ASF expectation values and write them to a CSV file.
    print *
    print '(a)', ">>> Evaluating ASF expectation values:"

    open(newunit=fh, file=trim(state)//".asfenergies.csv", action='write')
    call write_csv_header(fh)

    do k = 1, NVEC
        call eval_asfs(hamiltonian, k, hij)

        ! Write to the CSV file
        call write_csv_asf_line(fh, hij)

        ! Also print them out
        print '(80("="))'
        print '(">>> ASF #",i0)', k
        print '(80("-"))'
        print '(a1,a10,a18)',    " ", "Type", "<H_i> (Ha)"
        print '(a1,a10,e18.8)', " ", "DC", hij%diracpot + hij%coulomb
        print '(a1,a10,e18.8)', " ", "Breit", hij%breit
        print '(a1,a10,e18.8)', ">", "DC+B", hij%diracpot + hij%coulomb + hij%breit
        print '(a1,a10,e18.8)', " ", "QED VP", hij%vp
        print '(a1,a10,e18.8)', " ", "NMS+SMS", hij%nms + hij%sms
        print '(a1,a10,e18.8)', " ", "SE(Mohr)", hij%se_mohr
        print '(a1,a10,e18.8)', ">", "Full", &
            hij%diracpot + hij%coulomb + hij%breit + hij%vp + hij%nms + hij%sms + hij%se_mohr
    enddo

    close(fh)

contains

    subroutine write_csv_header(unit)
        integer, intent(in) :: unit
        write(unit, '(a9)', advance='no') "asf_index"
        write(unit, '(3(",",a16))', advance='no') "diracpot", "coulomb", "breit"
        write(unit, '(1(",",a16))', advance='no') "vp"
        write(unit, '(2(",",a16))', advance='no') "nms", "sms"
        write(unit, '(1(",",a16))', advance='no') "sum_dcb"
        write(unit, '()')
    end subroutine write_csv_header

    subroutine write_csv_asf_line(unit, me)
        integer, intent(in) :: unit
        type(matrixelement), intent(in) :: me

        write(unit, '(i9)', advance='no') k
        write(unit, '(3(",",e16.8))', advance='no') me%diracpot, me%coulomb, me%breit
        write(unit, '(1(",",e16.8))', advance='no') me%vp
        write(unit, '(2(",",e16.8))', advance='no') me%nms, me%sms
        write(unit, '(1(",",e16.8))', advance='no') &
            me%diracpot + me%coulomb + me%breit ! sum_dcb
        write(unit, '()') ! write the newline
    end subroutine write_csv_asf_line

    subroutine check_file(filename)
        use grasp_system, only: file_exists
        character(*), intent(in) :: filename

        if(.not.file_exists(filename)) then
            print '("Missing file ",a)', filename
            error stop
        endif
    end subroutine check_file

    subroutine eval_asfs(hamiltonian, k, hk)
        use grasp_kinds, only: real64, dp
        use orb_C
        use prnt_C
        use syma_C
        use eigv_C
        use eigvec1_C
        implicit none

        type(matrixelement), dimension(:,:), intent(in) :: hamiltonian
        integer, intent(in) :: k
        type(matrixelement), intent(out) :: hk

        integer :: i, j
        real(real64) :: cicj

        hk%diracpot = 0.0_dp
        hk%coulomb = 0.0_dp
        hk%breit = 0.0_dp
        hk%vp = 0.0_dp
        hk%nms = 0.0_dp
        hk%sms = 0.0_dp
        hk%se_mohr = 0.0_dp
        do i = 1, NCF
            do j = 1, NCF
                cicj = asf_coefficient(k, i) * asf_coefficient(k, j)
                hk%diracpot = hk%diracpot + hamiltonian(i, j)%diracpot * cicj
                hk%coulomb = hk%coulomb + hamiltonian(i, j)%coulomb * cicj
                hk%breit = hk%breit + hamiltonian(i, j)%breit * cicj
                hk%vp  = hk%vp  + hamiltonian(i, j)%vp * cicj
                hk%nms = hk%nms + hamiltonian(i, j)%nms * cicj
                hk%sms = hk%sms + hamiltonian(i, j)%sms * cicj
                hk%se_mohr = hk%se_mohr + hamiltonian(i, j)%se_mohr * cicj
            enddo
        enddo
    end subroutine eval_asfs

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
        use grasp_kinds, only: real64, dp
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
        use grasp_kinds, only: real64

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
        use grasp_kinds, only: real64
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

end program matrixelements
