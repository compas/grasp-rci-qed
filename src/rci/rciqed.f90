!> Types and routines for RCI.
module grasp_rciqed
    implicit none

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

end module grasp_rciqed
