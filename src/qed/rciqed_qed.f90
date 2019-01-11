module grasp_rciqed_qed
    use grasp_rciqed_kinds, only: real64
    use parameter_def, only: NNNW
    implicit none

    !> Stores the self-energy shifts for orbitals.
    real(real64) :: slfint(NNNW)

contains

    subroutine init_vacuum_polarization
        use decide_C, only: LVP
        use grid_C, only: N, RP
        use ncdist_C, only: ZDIST
        use tatb_C, only: TB
        use vpilst_C, only: NVPI, FRSTVP
        use ncharg_I
        use vacpol_I

        LVP = .TRUE.

        ! From AUXBLK
        CALL NCHARG
        CALL VACPOL
        ZDIST(2:N) = TB(2:N)*RP(2:N)
        FRSTVP = .TRUE.
        NVPI = 0
    end subroutine init_vacuum_polarization

end module grasp_rciqed_qed
