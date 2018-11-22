!> Routines related to the calculation of the normal and specific mass shifts.
module grasp_rciqed_mass_shifts
    implicit none

contains

    !> Initializes the nuclear mass shift global state.
    !!
    !! If the nuclear mass (`def_C: EMN`) is set to `0.0`, no initialization
    !! will occur and `LNMS` and `LSMS` are set to `.FALSE.`.
    subroutine init_mass_shifts
        use decide_C, only: LNMS, LSMS
        use def_C, only: EMN
        use keilst_C, only: NKEI, FRSTKI
        use vinlst_C, only: NVINTI, FRSTVI

        ! From AUXBLK
        IF (EMN > 0.D0) THEN
            !ATWINV = 1.D0/EMN
            LNMS = .TRUE.
            FRSTKI = .TRUE.
            NKEI = 0

            LSMS = .TRUE.
            FRSTVI = .TRUE.
            NVINTI = 0
        ELSE
            ! atwinv will not be used
            LNMS = .FALSE.
            LSMS = .FALSE.
        ENDIF
    end subroutine init_mass_shifts

end module grasp_rciqed_mass_shifts
