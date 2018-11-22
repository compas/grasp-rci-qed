!> Routines related to calculating the matrix elements of the Breit operator.
module grasp_rciqed_breit
    implicit none

contains

    !> Initialize Breit-related global state.
    !!
    !! @param j2max This is the value determined by genintrk(), e.g. in init_rkintc()
    subroutine init_breit(j2max)
        use parameter_def, only: NNNW
        use bcore_C, only: ICORE
        use bilst_C, only: FIRST, NTPI
        use decide_C, only: LTRANS
        use iounit_c, only: ISTDE
        use orb_C, only: NW, NCF
        use wfac_C, only: WFACT
        use genintbreit1_I
        use genintbreit2_I
        use ichop_I

        integer, intent(in) :: j2max

        ! These are the MPI parameters that need to be passed to different
        ! routines. We use the single core values.
        integer, parameter :: myid = 0, nprocs = 1

        integer :: i, j, N

        ! We'll enable all parts of the Hamiltonian. E.g. AUXBLK relies on these
        ! flags to determine if certain things get initialized.
        LTRANS = .TRUE.

        ! WFACT is the Breit scale factor. We set it to the default 1e-6
        WFACT = 1d-6

        ! From rci3mpi_LINUX.f
        call genintbreit1(myid, nprocs, N, j2max)
        call genintbreit2(myid, nprocs, N, j2max)

        ! The Breit part of auxblk.f90 (IF (LTRANS) THEN ...)
        ! Check the maximum numbers of orbtitals allowed in brint.f
        I = NNNW
        SELECT CASE (J2MAX)
        CASE (11)
            I = 114
        CASE (12)
            I = 112
        CASE (13)
            I = 110
        CASE (14)
            I = 108
        CASE (15)
            I = 106
        CASE (16)
            I = 105
        CASE (17)
            I = 103
        CASE (18)
            I = 101
        CASE (19)
            I = 100
        CASE (21:)
            I = 90
        END SELECT

        IF (I < NW) THEN
            WRITE (ISTDE, *) 'In setham. The number of orbitals is too'
            WRITE (ISTDE, *) 'large for the brint routine'
            STOP
        ENDIF

        FIRST = .TRUE.
        NTPI = 0

        ! BREID uses COMMON/BCORE/ICORE. This initialization was in setham_gg.f
        outer: do i = 1, NW
            ICORE(i) = 0
            do j = 1, NCF
                ! ICHOP is a lib92 routine
                if(ICHOP(i,j) <= 0) cycle outer
            enddo
            ICORE(i) = 1
        enddo outer

    end subroutine init_breit

end module grasp_rciqed_breit
