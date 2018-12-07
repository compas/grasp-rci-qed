!> Types and routines for RCI.
module grasp_rciqed
    implicit none

    !> I/O unit of the `.res` file.
    !!
    !! The default value is `-1`, so that it would not be possible to confuse
    !! this with any other open unit (negative unit number are generally not
    !! allowed and `newunit=` never returns a `-1`; i.e. it is never possible
    !! to write to `unit = -1`).
    integer :: res_unit = -1

contains

    !> Initializes the global state for the RK integrals, i.e. DC matrix elements.
    subroutine init_rkintc(j2max)
        use coeils_C, only: NCOEI, FRSTCO
        use genintrk_I

        ! These are the MPI parameters that need to be passed to different
        ! routines. We use the single core values.
        integer, parameter :: myid = 0, nprocs = 1

        integer, intent(inout) :: j2max

        integer :: N

        ! AUXBLK needs the j2max variable to be set. This is apparently set by
        ! GENINTRKwrap (based on rci3mpi.f). But GENINTRKwrap apparently only
        ! "wraps" the real call to genintrk, and then distributes it for MPI.
        ! We don't need it here, so we can just call genintrk directly..
        !
        ! genintrk() has another output though -- N, the "number of integrals".
        ! Whatever the significance of that is...
        !
        ! The values for the first two arguments (myid, nprocs) come from the
        ! non-mpi RCI rci92.f: we're the head node (myid = 0) of just 1 process
        ! (nprocs = 1).
        CALL genintrk(myid, nprocs, N, j2max)

        ! From AUXBLK
        FRSTCO = .TRUE.
        NCOEI = 0

    end subroutine init_rkintc

    !> Writes a $(jobname).settings.toml file with the nuclear, grid and
    !! Hamiltonian settings that were used for the RCI run.
    subroutine write_job_toml(jobname)
        use decide_C, only: LTRANS, LNMS, LSMS, LVP, LSE
        use def_C, only: Z, EMN, AUMAMU, FMTOAU
        use grid_C, only: RNT, H, N
        use npar_C, only: NPARM, PARM

        character(len=*), intent(in) :: jobname
        character(len=:), allocatable :: tomlfile
        integer :: toml_unit

        tomlfile = trim(jobname) // ".settings.toml"

        open(newunit=toml_unit, file=tomlfile)
        write(toml_unit, '(a)') "[nucleus]"
        write(toml_unit, '(a,f0.16)') "  Z = ", Z
        write(toml_unit, '(a,f0.16)') "  atomic_mass_amu = ", EMN * AUMAMU
        if(NPARM == 0) then
            write(toml_unit, '(a)') "  nuclear_model = ""point"""
        elseif(NPARM == 2) then
            write(toml_unit, '(a)') "  nuclear_model = ""fermi"""
            write(toml_unit, '(a,e22.16)') "  fermi_a = ", PARM(2) / FMTOAU
            write(toml_unit, '(a,e22.16)') "  fermi_c = ", PARM(1) / FMTOAU
        else
            print *, "ERROR: Bad nuclear model; NPARM = ", NPARM
            write(toml_unit, '(a,i0)') "  nuclear_model = ", NPARM
        endif

        write(toml_unit, '(a)') "[grid]"
        write(toml_unit, '(a,e22.16,a)') "  R = ", RNT, " # = 2e-6 / Z"
        write(toml_unit, '(a,e22.16)') "  H = ", H
        write(toml_unit, '(a,i0)') "  N = ", N

        write(toml_unit, '(a)') "[hamiltonian]"
        write(toml_unit, '(a)') "  # Contains the following booleans:"
        write(toml_unit, '(a)') "  #   breit, nms, sms, qed_vp, qed_se"
        write(toml_unit, '(a)') "  # Indicating if the corresponding part of the Hamiltonian"
        write(toml_unit, '(a)') "  # was enabled. If it is omitted, it is assumed to have been off."
        if(LTRANS) write(toml_unit, '(a)') "  breit = true"
        if(LNMS) write(toml_unit, '(a)') "  nms = true"
        if(LSMS) write(toml_unit, '(a)') "  sms = true"
        if(LVP) write(toml_unit, '(a)') "  qed_vp = true"
        if(LSE) write(toml_unit, '(a)') "  qed_se = true"

        close(toml_unit)

        print *, tomlfile
    end subroutine write_job_toml

end module grasp_rciqed
