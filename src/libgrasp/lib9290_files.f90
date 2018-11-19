!> Routines for initializing GRASP global state from input files (`isodata`,
!! CSLs, wavefunction files, mixing files).
module grasp_lib9290_files
    implicit none

contains
    ! NOTE -- the following applies to all the routines:
    !
    ! These are the lib92 initialization calls. Not really sure which
    ! ones are necessary and how they depend on each other. Getting this
    ! to work was a trial and error matter.
    !
    ! But the basic ordering was taken from rci's rci92.f main file.
    !   - SETDBG was ignored, since it's a RCI-specific routine and did
    !     not look relevant.
    !   - SETSUM was ignored as well.. "Open summary file" does not sound
    !     relevant either.

    !> Loads a configuration state list (CSL) files (typically either `rcsf.inp`
    !! or `*.c`).
    !!
    !! TODO: Document prerequisites (what needs to be loaded before this can be
    !! called?).
    subroutine load_csl(cfile)
        use grasp_lib9290_csls, only: count_blocks
        use memory_man
        use alcbuf_I
        use setcsll_I
        use lodcsh_I
        use lodcsh2_I
        use setqic_I
        use radgrd_I
        use parameter_def, only: NNNW
        use hblock_C, only: NBLOCK, NCFBLK
        use orb_C, only: NCF, IQA
        use stat_C, only: JQSA, JCUPA

        character(*), intent(in) ::  cfile
        character(:), allocatable :: cfile_local
        integer :: ncore, nblocks
        integer :: status

        character(8), allocatable :: idblk(:)

        ! The following is based on the SETCSL routine from RCI.
        call count_blocks(cfile, nblocks, status)
        allocate(idblk(nblocks))
        !CALL ALLOC(pncfblk, nblocks, 4)
        CALL ALLOC(ncfblk, nblocks, 'ncfblk', 'lib92_init_cw')
        ! TODO: SETCSLL is pretty trivial. It could be re-written to be more sane.
        ! NOTE: SETCSL always allocating N+1 elements.. I decided to go with N,
        ! but unsure if it has some reason maybe?
        ! NOTE: In SETCSL the second to last argument (NCF) was called ncftot in
        ! the ORB2 common block.
        cfile_local = cfile
        CALL SETCSLL(21, cfile_local, nblocks, NBLOCK, NCFBLK, NCF, idblk)
        ! As we already have `nblock` elements in `ncfblk`, there is no need
        ! to re-allocate, as was done in SETCSL.
        REWIND (21)
        READ (21,*) ! WTH? But LODCSH fails without this.
        CALL LODCSH (21, ncore)
        ! Ok, LODCSH only loads the _header_ of the CSL file. Among probably other
        ! stuff, it leaves the JCUPA variable undefined, causing segfaults.
        !
        ! Ok.. LODCSH2 seems like it might load the CSFs in. It requires
        !   - ORB2/PNTRIQ
        !   - STAT/PNTJQS
        !   - STAT/PNJCUP
        ! to be allocated though. The allocation was taken from matrix.f.
        !
        ! NOTE: to get it working, before I figure LODCSH2, I also tried
        !   call LODCSL(ncore)
        !   call LODCSL(21, ncore)
        !   call LODCSL2(21, ncore, 1)
        !
        ! These are from the GRASP2K implementation..
        ! CALL ALLOC(PNTRIQ, NNNWP*NCF, 4)
        ! CALL ALLOC(PNTJQS, NNNWP*NCF*3, 4)
        ! CALL ALLOC(PNJCUP, NNNWP*NCF, 4)
        !
        ! There are from matrix.f
        CALL ALLOC(IQA, NNNW, ncf, 'IQA', 'MATRIX')
        CALL ALLOC(JQSA, NNNW, 3, ncf, 'JQSA', 'MATRIX')
        CALL ALLOC(JCUPA, NNNW, ncf, 'JCUPA', 'MATRIX')
        ! CALL ALLOC (SLF_EN,ncf,'SLF_EN', 'MATRIX')
        ! CALL ALLOC (UCF, nw,'UCF', 'MATRIX')
        CALL LODCSH2(21, ncore, -119)
    end subroutine load_csl

    !> Loads the nuclear data and initializes the grid and the nuclear model.
    subroutine load_isodata(isofile)
        use grasp_lib9290_init, only: lib9290_init_grid
        use def_C, only: CVAC, C, Z
        use nucpot_I
        use setiso_I

        character(*), intent(in) :: isofile

        ! What follows is from GETCID (RCI92 -> SETRES -> GETCID)
        call SETISO(isofile)
        C = CVAC
        call lib9290_init_grid(Z)
        CALL NUCPOT
    end subroutine load_isodata

    !> Loads the radial orbitals.
    subroutine load_orbitals(wfile)
        use setrwfa_I

        character(*), intent(in) :: wfile

        ! What follows is from GETCID (RCI92 -> SETRES -> GETCID)
        call SETRWFA(wfile)
    end subroutine load_orbitals

    !> Loads a mixing files (typically `*.m` or `*.cm`).
    subroutine load_mixing(sumfile)
        use getmixblock_I
        character(*), intent(in) :: sumfile
        call getmixblock(sumfile) ! TODO: implicit local routine
    end subroutine load_mixing

end module grasp_lib9290_files
