!> Private module containing GRASP `lib92` common block definitions used by the
!! routines in `g2k_lib92`.
!!
!! A separate module is used so that it would be possible to use `implicit none`
!! in the main function body.
!!
!! The module __should not__ be used outside of the `g2k_lib92.f90` file.
!!
! module g2k_lib92_common
!     use g2k_parameters
!     implicit real*8 (a-h, o-z)
!
!     private
!
!     pointer (PNCFBLK,NCFBLK(1)) ! originally: POINTER (pncfblk, ncfblk(0:*))
!
!     COMMON/ORB2/NCF,NW,PNTRIQ &
!           /DEF2/C &
!           /DEF4/ACCY,NSCF,NSIC,NSOLV &
!           /DEF9/CVAC,PI &
!           /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N &
!           /HBLOCK/NBLOCK,PNCFBLK &
!           /STAT/PNTJQS,PNJCUP
!
!     public :: NNNWP
!     public :: C, CVAC, ACCY, H
!     public :: NCF, PNTRIQ
!     public :: PNCFBLK, NBLOCK, NCFBLK
!     public :: PNTJQS, PNJCUP
! end module g2k_lib92_common


!> Wrapper methods used to interact with the `lib92` library, collected here for
!! potential future reuse.
module g2k_lib92
    use orb_C
    use def_C
    use grid_C
    use hblock_C
    use stat_C
    implicit none

contains

    !> Initializes the GRASP common blocks based on the input files.
    !!
    !! @param isofile The isotope and grid data file, conventionally called `isodata`.
    !! @param cfile The configuration state list (CSL) file (usually `*.c`).
    !! @param wfile Orbital wavefunction file (usually `*.w`).
    !!
    subroutine lib92_init_cw(isofile, cfile, wfile)
        use g2k_csls, only: count_blocks
        use grasp_lib9290_init
        !use g2k_lib92_common
        use memory_man
        use alcbuf_I
        use setcsll_I
        use lodcsh_I
        use lodcsh2_I
        use setiso_I
        use setqic_I
        use radgrd_I
        use nucpot_I
        use setrwfa_I
        use factt_I
        use setmc_I
        use setcon_I

        character(*), intent(in) :: isofile, cfile, wfile
        character(:), allocatable :: cfile_local
        integer :: ncore, nblocks
        integer :: status

        character(8), allocatable :: idblk(:)

        ! These are the lib92 initialization calls. Not really sure which
        ! ones are necessary and how they depend on each other. Getting this
        ! to work was a trial and error matter.
        !
        ! But the basic ordering was taken from rci's rci92.f main file.
        !   - SETDBG was ignored, since it's a RCI-specific routine and did
        !     not look relevant.
        !   - SETSUM was ignored as well.. "Open summary file" does not sound
        !     relevant either.
        call SETMC ! Perform machine- and installation-dependent setup
        call SETCON ! Set up the physical constants

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

        ! What follows is from GETCID (RCI92 -> SETRES -> GETCID)
        call SETISO(isofile)
        C = CVAC
        call lib9290_init_grid(Z)
        CALL NUCPOT
        call SETRWFA(wfile)

        call lib9290_init_rkco_gg
    end

    subroutine lib92_init_csls(cfile)
        use g2k_csls, only: count_blocks
        !use g2k_lib92_common
        use memory_man
        use alcbuf_I
        use setcsll_I
        use lodcsh_I
        use lodcsh2_I
        use setiso_I
        use setqic_I
        use radgrd_I
        use nucpot_I
        use setrwfa_I
        use factt_I
        use setmc_I
        use setcon_I

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
    end subroutine lib92_init_csls

    subroutine lib92_init_mixing(sumfile)
        use getmixblock_I
        character(*), intent(in) :: sumfile
        call getmixblock(sumfile) ! TODO: implicit local routine
    end
end module g2k_lib92
