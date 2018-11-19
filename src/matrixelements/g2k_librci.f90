!> Private module containing GRASP common block definitions used by the
!! `g2k_lib92::rcicommon_*` routines.
!!
!! A separate module is used so that it would be possible to use `implicit none`
!! in the main function body.
!!
!! The module __should not__ be used outside of the `g2k_lib92.f90` file.
!!
! module g2k_librci_common
!     use g2k_parameters
!     implicit real*8 (a-h, o-z)
!
!     private
!
!     ! ICHOP is from lib92
!     integer, external :: ICHOP
!
!     LOGICAL LFORDR,LTRANS,LVP,LSE,LNMS,LSMS
!     COMMON/DECIDE/LFORDR,LTRANS,LVP,LSE,LNMS,LSMS &
!           /BCORE/ICORE(NNNW) &
!           /ORB2/NCF,NW,PNTRIQ &
!           /WFAC/WFACT
!
!     public :: LTRANS, LVP, LSE, LNMS, LSMS
!     public :: ICORE, ICHOP
!     public :: NW, NCF
!     public :: WFACT
! end module g2k_librci_common

!> Wrapper methods used to interact with the `lib92` library, collected here for
!! potential future reuse.
module g2k_librci
    use grasp_kinds, only: real64
    !use g2k_librci_common, only: NW
    use orb_C, only: NW
    implicit none

    type hamiltonian_cache
        real(real64), dimension(:,:,:), allocatable :: sematrices
    end type hamiltonian_cache

    type matrixelement
        real(real64) :: diracpot, coulomb, breit, &
            vp, se(4), &
            nms, sms
    end type matrixelement

contains

    !> ...
    subroutine allocate_sematrices(cache)
        !use g2k_librci_common, only: NW
        type(hamiltonian_cache), intent(inout) :: cache
        integer :: setype

        allocate(cache%sematrices(4, NW, NW))
        do setype = 1, 4
            !call qedse(setype - 1, cache%sematrices(setype, :, :))
        enddo
    end subroutine allocate_sematrices

    !> Initializes the global state for the RK integrals, i.e. DC matrix elements.
    subroutine init_rkintc(j2max)
        !use g2k_parameters, only: real64
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
        print '("N = ",i0)', N
        print '("j2max = ",i0)', j2max

        ! From AUXBLK
        FRSTCO = .TRUE.
        NCOEI = 0

    end subroutine init_rkintc

end module g2k_librci
