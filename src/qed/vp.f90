module grasp_rciqed_qed_vp
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNNP
    implicit none

    real(real64), dimension(:), allocatable :: vac2p4

contains

    !> Populates `matrix` with the QED vacuum polarization matrix elements in the orbital
    !! basis.
    !!
    !! The matrix element values are integrated using the `vac2p4` global array from this
    !! module.
    !!
    !! @param matrix An `NW x NW` `real64` array for storing the matrix elements.
    subroutine qedvp(matrix)
        use orb_C, only: NW, NAK
        ! Arguments:
        real(real64), intent(out) :: matrix(:, :)
        ! Local variables:
        real(real64) :: vpij
        integer :: k, l

        do k = 1, NW
            do l = k, NW
                if(NAK(k) == NAK(l)) then
                    vpij = potential(vac2p4, k, l)
                else
                    vpij = 0.0_dp
                endif
                matrix(k, l) = vpij
                matrix(l, k) = vpij
            enddo
        enddo
    end

    !> Integrates the generic potential, stored in the array `v`, with orbitals `k1` and
    !! `k2` to obtain a single particle matrix element correponding to those orbitals.
    !!
    !! Internally, it uses the `QUAD` routine to perform the integration on the standard
    !! GRASP grid. It will multiply the value of the potential with `RP`, so the user should
    !! do that themselves (i.e. internally `TA ~ V * RP`).
    !!
    !! @param v An array containing the potential, represented on the GRASP grid.
    !! @param k1,k2 Indices of the orbitals.
    !!
    !! @return The value of the single particle matrix element for orbitals `k1` and `k2`.
    function potential(v, k1, k2)
        use grid_C, only: N, RP
        use tatb_C, only: MTP, TA
        use wave_C, only: PF, QF
        use quad_I
        ! Arguments:
        real(real64), dimension(:), intent(in) :: v
        integer, intent(in) :: k1, k2
        real(real64) :: potential
        ! Local variables:
        integer :: i

        do i = 1, N
            TA(i) = (PF(i,k1)*PF(i,k2) + QF(i,k1)*QF(i,k2)) * RP(i) * v(i)
        end do
        MTP = N
        call QUAD(potential)
    end

    ! Legacy GRASP92 routines. Modifications to these are minimal, primarily just to make
    ! them independent of old GRASP common blocks / _c modules.
    include 'legacy/funk.f90'
    include 'legacy/funl.f90'

end module grasp_rciqed_qed_vp
