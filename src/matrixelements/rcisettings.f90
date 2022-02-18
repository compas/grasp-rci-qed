module grasp_rciqed_rcisettings
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

    public :: rcisettings, read_settings_toml, write_settings_toml

    type rcisettings
        ! Nuclear parameters
        real(real64) :: Z, atomic_mass_amu
        !NPARM, PARM(2)
        ! Grid parameters
        !real(real64) :: RNT, H, N
        ! Hamiltonian parameters
        character(len=:), allocatable :: breit_mode
        logical :: nms_enabled, sms_enabled, qed_vp_enabled
        ! For self-energy, we store:
        integer :: qed_se, qed_se_hydrogenic_cutoff
    end type rcisettings

contains

    !> Read a `.setting.toml` for an RCI job into a `rcisettings` object.
    !!
    !! @param jobname Name of the RCI job. File assumed to be called `$(jobname).settings.toml`.
    !! @param settings `rcisettings` instance to be populated by the values.
    !! @return Returns `.true.` it the parsing was successful, and `.false.` if
    !!         there was a problem.
    function read_settings_toml(jobname, settings)
        use grasp_rciqed_system, only: file_exists
        use grasp_rciqed_toml
        use grasp_rciqed_qed, only: nsetypes, setypes_short

        character(len=*), intent(in) :: jobname
        type(rcisettings), intent(out) :: settings
        logical :: read_settings_toml

        character(len=:), allocatable :: tomlfile, setype_string
        type(cpptoml_table) :: table
        integer :: i

        tomlfile = trim(jobname) // '.settings.toml'
        if(.not.file_exists(tomlfile)) then
            print '("WARNING[read_settings_toml]: ",a," does not exist.")', tomlfile
            read_settings_toml = .false.
            return
        endif
        call parse_file(tomlfile, table)

        if(.not.get_double(table, "nucleus.Z", settings%Z)) then
            print '("WARNING[read_settings_toml]: nucleus.Z missing.")'
            read_settings_toml = .false.
            return
        endif

        if(.not.get_double(table, "nucleus.atomic_mass_amu", settings%atomic_mass_amu)) then
            print '("WARNING[read_settings_toml]: nucleus.atomic_mass_amu missing.")'
            read_settings_toml = .false.
            return
        endif

        ! TODO: Load nuclear model and GRID too!

        call get_string_default(table, "hamiltonian.breit", "n", settings%breit_mode)
        call get_logical_default(table, "hamiltonian.nms", .false., settings%nms_enabled)
        call get_logical_default(table, "hamiltonian.sms", .false., settings%sms_enabled)
        call get_logical_default(table, "hamiltonian.qed_vp", .false., settings%qed_vp_enabled)

        ! Parse the self-energy
        settings%qed_se = 0 ! If hamiltonian.qed_se is missing => QED SE was disabled
        if(get_string(table, "hamiltonian.qed_se", setype_string)) then
            do i = 1, nsetypes
                if(trim(setypes_short(i)) == trim(setype_string)) then
                    settings%qed_se = i
                    exit
                endif
            enddo
            if(settings%qed_se == 0) then
                print '("WARNING[read_settings_toml]: hamiltonian.qed_se has bad value.")'
                print '("  Contains: ",a)', trim(setype_string)
                print '("  Expected one of: ",10(a,:", "))', (trim(setypes_short(i)), i = 1,4)
                read_settings_toml = .false.
                return
            endif
        endif
        call get_integer_default(table, "hamiltonian.qed_se_hydrogenic_cutoff", -1, settings%qed_se_hydrogenic_cutoff)

        ! If we didn't return by now, it all went well
        read_settings_toml = .true.
    end function read_settings_toml

    !> Writes a $(jobname).settings.toml file with the nuclear, grid and
    !! Hamiltonian settings that were used for the RCI run.
    subroutine write_settings_toml(jobname, breit_mode, setype)
        use decide_C, only: LTRANS, LNMS, LSMS, LVP, LSE
        use def_C, only: Z, EMN, AUMAMU, FMTOAU
        use grid_C, only: RNT, H, N
        use npar_C, only: NPARM, PARM
        use qedcut_C, only: NQEDCUT, NQEDMAX
        use grasp_rciqed_qed, only: nsetypes, setypes_short

        character(len=*), intent(in) :: jobname
        ! TODO: passing setype and breit_mopde explicitly here like this, even though
        ! otherwise we rely on the global common block modules is a bit of an inconsistent
        ! hack to resolve a circular dependency problem. This should be fixed at some point.
        character(len=*), intent(in) :: breit_mode
        integer, intent(in) :: setype

        character(len=:), allocatable :: tomlfile
        integer :: toml_unit

        tomlfile = trim(jobname) // ".settings.toml"

        open(newunit=toml_unit, file=tomlfile)
        write(toml_unit, '(a)') "[nucleus]"
        write(toml_unit, '(a,f0.16)') "  Z = ", Z
        if(EMN == 0) then
            ! Fortran prints '.000000000' if the floating point number is exactly zero,
            ! which breaks the TOML file (not a valid float literal). So we need to
            ! special-case that.
            write(toml_unit, '(a)') "  atomic_mass_amu = 0.0"
        else
            write(toml_unit, '(a,f0.16)') "  atomic_mass_amu = ", EMN * AUMAMU
        endif
        if(NPARM == 0) then
            write(toml_unit, '(a)') "  nuclear_model = ""point"""
        elseif(NPARM == 2) then
            write(toml_unit, '(a)') "  nuclear_model = ""fermi"""
            write(toml_unit, '(a,e22.16)') "  #fermi_a: ", PARM(2) / FMTOAU
            call write_toml_expfloat(toml_unit, "fermi_a", PARM(2) / FMTOAU)
            write(toml_unit, '(a,e22.16)') "  #fermi_c: ", PARM(1) / FMTOAU
            call write_toml_expfloat(toml_unit, "fermi_c", PARM(1) / FMTOAU)
        else
            print *, "WARNING: Bad nuclear model; NPARM = ", NPARM
            write(toml_unit, '(a,i0)') "  nuclear_model = ", NPARM
        endif

        write(toml_unit, '(a)') "[grid]"
        call write_toml_expfloat(toml_unit, "RNT", RNT)
        write(toml_unit, '(a,e22.16,a)') "  #R: ", RNT, " # = 2e-6 / Z"
        call write_toml_expfloat(toml_unit, "H", H)
        write(toml_unit, '(a,e22.16)') "  #H: ", H
        write(toml_unit, '(a,i0)') "  N = ", N

        write(toml_unit, '(a)') "[hamiltonian]"
        write(toml_unit, '(a)') "  # Contains the following booleans:"
        write(toml_unit, '(a)') "  #   qed_vp, nms, sms"
        write(toml_unit, '(a)') "  # Each indicates if the corresponding part of the Hamiltonian"
        write(toml_unit, '(a)') "  # was enabled. If it is omitted, it is assumed to have been off."
        write(toml_unit, '(a)') "  #"
        write(toml_unit, '(a)') "  # Can also contain breit (string) which indicates the mode of the"
        write(toml_unit, '(a)') "  # transverse photon operator. Possible values are 'breit', 'specorbs'"
        write(toml_unit, '(a)') "  # and 'full'. If set to 'n', the operator is completely disabled."
        write(toml_unit, '(a)') "  #"
        write(toml_unit, '(a)') "  # May also contain qed_se (string), which indicates that a particular"
        write(toml_unit, '(a)') "  # QED self-energy operator was also included in the Hamiltonian."
        write(toml_unit, '(a)') "  # The possible values are: 'hydrogenic', 'qedmod', 'flambaum' and 'pyykkoe'"
        write(toml_unit, '(a)') "  #"
        write(toml_unit, '(a)') "  # Finally, for the hydrogenic QED self-energy, qed_se_hydrogenic_cutoff"
        write(toml_unit, '(a)') "  # (integer) may be defined, which sets the n-quantum number cutoff. If not"
        write(toml_unit, '(a)') "  # present, the implementation default is used."
        if(LTRANS) then
            write(toml_unit, '("  breit = """,a,"""")') trim(breit_mode)
        else
            write(toml_unit, '("  breit = ""n""")')
        endif
        if(LNMS) write(toml_unit, '(a)') "  nms = true"
        if(LSMS) write(toml_unit, '(a)') "  sms = true"
        if(LVP) write(toml_unit, '(a)') "  qed_vp = true"
        if(LSE) then
            if(setype < 1 .or. setype > nsetypes) then
                print *, "ERROR: Bad self-energy type."
                error stop
            endif
            write(toml_unit, '("  qed_se = """,a,"""")') trim(setypes_short(setype))
            if(setype == 1 .and. NQEDCUT == 1) then
                write(toml_unit, '(a,i0)') "  qed_se_hydrogenic_cutoff = ", NQEDMAX
            endif
        endif

        close(toml_unit)

        print *, tomlfile
    end subroutine write_settings_toml

    !> Writes a TOML key-value line for a floating point number in exponent notation.
    !!
    !! The problem is that with the `E` Fortran edit descriptors you end up with
    !! a leading zero in the exponent, which is not allowed by the TOML format
    !! specification. To the best of my knowledge, you can not fix this with
    !! format descriptions in a generic way (you _can_ set the width of the
    !! exponent to 1 character, but then writing numbers with exponents
    !! \f$ |\textrm{exp}| \geq 10 \f$ will fail).
    subroutine write_toml_expfloat(unit, label, value)

        integer, intent(in) :: unit
        character(*), intent(in) :: label
        real(real64), intent(in) :: value

        logical :: sign ! .true. => positive
        integer :: exponent
        real(real64) :: significand

        if(value == 0.0_dp) then
            sign = .true.
            exponent = 0
            significand = 0.0_dp
        else
            sign = (value >= 0.0_dp)
            exponent = 1 + floor(log10(abs(value)))
            significand = abs(value) / (10.0_dp ** exponent)
        endif

        if(sign) then
            write(unit, '("  ",a," = 0",f0.16,"e",i0)') label, significand, exponent
        else
            write(unit, '("  ",a," = -0",f0.16,"e",i0)') label, significand, exponent
        endif
    end subroutine write_toml_expfloat

end module grasp_rciqed_rcisettings
