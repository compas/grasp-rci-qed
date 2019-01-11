!> Wrapper methods used to interact with the `lib92` library, collected here for
!! potential future reuse.
module grasp_lib9290
    implicit none

contains

    !> Initializes the GRASP global state based on the input files.
    !!
    !! This attempts to be the one-stop-shop for `lib9290` initialization.
    !!
    !! @param isofile The isotope and grid data file, conventionally called `isodata`.
    !! @param cfile The configuration state list (CSL) file (usually `*.c`).
    !! @param wfile Orbital wavefunction file (usually `*.w`).
    !!
    subroutine init_isocw_full(isofile, cfile, wfile)
        use grasp_rciqed_lib9290_init
        use grasp_rciqed_lib9290_files

        character(*), intent(in) :: isofile, cfile, wfile

        call lib9290_init_constants ! physical and machine constants
        call load_csl(cfile)
        call load_isodata(isofile)
        call load_orbitals(wfile)
        call lib9290_init_rkco_gg ! TODO: does this need to be here?
    end

end module grasp_lib9290
