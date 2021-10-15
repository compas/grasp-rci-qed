module grasptest_qedvp_legacy
    implicit none

contains

    !> This routine sets up the necessary global state for the legacy QED VP routines.
    !! It should be called after qedvp_init.
    subroutine qedvp_legacy_init
        use grasp_rciqed_qed_vp, only: qedvp_initialized
        use grasptest_vpilst_C, only: NVPI, FRSTVP

        if(.not.qedvp_initialized) then
            print *, "ERROR(grasp_rciqed_qed_vp): GRASP VP global state not initialized."
            print *, "qedvp_init() has to be called before the other routines can be used."
            stop 1
        end if

        ! The following two variables are used by the VPINT legacy routine, to manage the
        ! dynamic storage of the matrix elements. These values were used to be set in AUXBLK.
        FRSTVP = .TRUE.
        NVPI = 0

    end

end
