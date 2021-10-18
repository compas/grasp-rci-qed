module grasp_rciqed_qed_vp
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNNP
    implicit none

    real(real64), allocatable :: vp_vac2(:), vp_vac4(:)
    logical :: qedvp_initialized  = .false.

    ! Legacy global array from the ncdist_C module:
    real(real64), DIMENSION(NNNP) :: ZDIST

contains

    !> Initializes the global state for the vacuum polarization calculations.
    !!
    !! If the global state is already initialized, the routine just prints an error message
    !! and does not do anything.
    subroutine qedvp_init
        use decide_C, only: LVP
        use grid_C, only: N, R, RP
        use tatb_C, only: MTP, TB
        use ncharg_I
        use vac2_I
        use vac4_I

        if(qedvp_initialized) then
            print *, "ERROR(grasp_rciqed_qed_vp/qedvp_init): GRASP VP global state already initialized."
            return
        end if

        LVP = .TRUE.

        ! From AUXBLK. NCHARG set up the nuclear charge distrbution in ZDIST. If we're dealing
        ! with the PNC, then the routine actually does nothing, ZDIST stays zero and is
        ! ignored in VAC2 and VAC4.
        CALL NCHARG
        ! VACPOL scaled ZDIST before calling VAC2 and VAC4, so we also need to do that here.
        ! MTP (from tatb_C) is set by NCHARG.
        ZDIST(:MTP) = ZDIST(:MTP)*R(:MTP)*RP(:MTP)
        ! The VP potentials are calculated in VAC2 and VAC4, which both write it into the
        ! TB array.
        call VAC2
        allocate(vp_vac2(N))
        vp_vac2(:N) = TB(:N)
        ! VAC4, however, does not overwrite TB, but rather adds to it (in the original
        ! implementation it was adding the VAC4 potential on top of VAC2). However, as we
        ! just want to get the VAC4 potential contribution, we set TB back to zero here
        ! first.
        TB(:N) = 0.0_dp
        call VAC4
        allocate(vp_vac4(N))
        vp_vac4(:N) = TB(:N)

        qedvp_initialized = .true.
    end

    !> Populates `matrix` with the QED vacuum polarization matrix elements in the orbital
    !! basis.
    !!
    !! The matrix element values are integrated using the `vac2p4` global array from this
    !! module, which contains both the Uehling and Källen-Sabry terms.
    !!
    !! @param matrix An `NW x NW` `real64` array for storing the matrix elements.
    subroutine qedvp(matrix)
        use orb_C, only: NW
        ! Arguments:
        real(real64), intent(out) :: matrix(:, :)
        ! Local variables:
        real(real64) :: vpij
        integer :: k, l

        if(.not.qedvp_initialized) then
            print *, "ERROR(grasp_rciqed_qed_vp): GRASP VP global state not initialized."
            print *, "qedvp_init() has to be called before the other routines can be used."
            error stop
        end if

        do k = 1, NW
            do l = k, NW
                vpij = qedvp_kl(k, l)
                matrix(k, l) = vpij
                matrix(l, k) = vpij
            enddo
        enddo
    end

    !> Returns the value of the QED VP matrix element for orbitals `k` and `l`.
    !!
    !! The matrix element values are integrated using the `vp_vac2` and `vp_vac4` global
    !! arrays from this module, which contain the Uehling and Källen-Sabry terms,
    !! respectively.
    !!
    !! @param k, l Indices of the orbital for which the VP integral is calculated for.
    !!
    !! Note: no error checking is done, neither for whether VP has been initialized with
    !! `qedvp_init`, nor whether the orbital indices are within bounds.
    function qedvp_kl(k, l)
        use orb_C, only: NAK
        ! Arguments:
        integer, intent(in) :: k, l
        real(real64) :: qedvp_kl

        ! Note: potential also does this check, but if we do it here, we can avoid
        ! evaluating the vp_vac2 + vp_vac4 sum on zero elements.
        if(NAK(k) == NAK(l)) then
            qedvp_kl = potential(vp_vac2 + vp_vac4, k, l)
        else
            qedvp_kl = 0.0_dp
        end if
    end

    !> Integrates the generic potential, stored in the array `v`, with orbitals `k1` and
    !! `k2` to obtain a single particle matrix element correponding to those orbitals.
    !!
    !! Internally, it uses the `QUAD` routine to perform the integration on the standard
    !! GRASP grid. It will multiply the value of the potential with `RP`, so the user should
    !! not do that themselves (i.e. internally `TA ~ V * RP`).
    !!
    !! @param v An array containing the potential, represented on the GRASP grid.
    !! @param k1,k2 Indices of the orbitals.
    !!
    !! @return The value of the single particle matrix element for orbitals `k1` and `k2`.
    function potential(v, k1, k2)
        use grid_C, only: N, RP
        use tatb_C, only: MTP, TA
        use orb_C, only: NAK
        use wave_C, only: PF, QF, MF
        use quad_I
        ! Arguments:
        real(real64), dimension(:), intent(in) :: v
        integer, intent(in) :: k1, k2
        real(real64) :: potential
        ! Local variables:
        integer :: i

        ! If the orbitals are from different angular momenta, then the matrix element of the
        ! potential is zero due to symmetry.
        if(NAK(k1) /= NAK(k2)) then
            potential = 0.0_dp
            return
        end if
        ! MF appears to store the MTP for the orbitals, and we'll use that to determine the
        ! MTP of the integral. This is consistent with the old VPINTF routine.
        MTP = min(MF(k1), MF(k2))
        do i = 1, MTP
            TA(i) = (PF(i,k1)*PF(i,k2) + QF(i,k1)*QF(i,k2)) * RP(i) * v(i)
        end do
        call QUAD(potential)
    end

    ! Legacy GRASP92 routines. Modifications to these are minimal, primarily just to make
    ! them independent of old GRASP common blocks / _c modules.
    include 'legacy/funk.f90'
    include 'legacy/funl.f90'

end module grasp_rciqed_qed_vp
