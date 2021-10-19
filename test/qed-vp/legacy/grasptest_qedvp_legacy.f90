module grasptest_qedvp_legacy
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNNP
    implicit none

    ! Legacy global array from the ncdist_C module. Note: there is an array with
    ! the same name in grasp_rciqed_qed_vp. This one here specifically is used
    ! in the legacy VPINTF routine purely for testing purposes.
    real(real64), DIMENSION(NNNP) :: ZDIST

contains

    !> This routine sets up the necessary global state for the legacy QED VP routines.
    !! It should be called after qedvp_init.
    subroutine qedvp_legacy_init
        use grid_C, only: N, RP
        use grasp_rciqed_qed_vp, only: qedvp_initialized, vp_potential
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

        ! VACPOL puts the vacuum polarization potentials into TB, but the VPINT(F) routines
        ! use ZDIST to actually evaluate the matrix elements. This scaling of the TB array
        ! used to be in AUXBLK.
        ZDIST(1) = 0.0_dp
        ZDIST(2:N) = (vp_potential(1, 2:N) + vp_potential(2, 2:N)) * RP(2:N)

    end

end
