program hydrogenic_vp
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use grasp_rciqed_qed_vp, only: qedvp_init, qedvp, qedvp_kl
    use grasptest_lib9290_setup
    use grasptest_testing
    use grasptest_utilities
    use grasptest_qedvp_legacy, only: qedvp_legacy_init
    use orb_C
    implicit none

    ! Tolerances:
    real(real64), parameter :: atol = 1e-15_dp, rtol = 1e-15_dp

    type orbital_reference
        integer :: orbital_n, orbital_kappa
        real(real64) :: vp, vp_tol ! self-energy and its tolerance
    end type orbital_reference

    ! This reference data was simply calculated using the same code it is testing.
    ! I.e. it is not from a separate reference or an independent calculation and
    ! therefore one should not assume that these numbers are necessarily correct.
    ! They are here to catch any change to the diagonal element values.
    ! The reference values are accurate up to 10^{-8}.
    type(orbital_reference), parameter, dimension(*) :: orbitals = (/ &
        orbital_reference(1, -1, -3.1703841900e-03, 1d-7), & ! 1s
        orbital_reference(2, -1, -4.0118488628e-04, 1d-7), & ! 2s
        orbital_reference(2,  1, -1.8052645141e-06, 1d-7), & ! 2p-
        orbital_reference(2, -2, -3.4397524133e-07, 1d-7), & ! 2p
        orbital_reference(3,  2, -1.1241603499e-10, 1d-7), & ! 3d-
        orbital_reference(3, -3, -3.3478594312e-11, 1d-7), & ! 3d
        orbital_reference(4,  3, -6.9328886357e-15, 1d-7), & ! 4f-
        orbital_reference(5, -4, -2.0898993567e-15, 1d-7)  & ! 5f
    /)
    type(orbital_reference) :: orb

    integer :: n, k, l
    real(real64) :: vp_value
    real(real64), dimension(:,:), allocatable :: vp_matrix
    character(len=:), allocatable :: label
    logical :: success = .true.

    call setup(18.0_dp, 39.9623831225_dp) ! Ar-40
    call allocate_hydrogenic_orbitals

    call qedvp_init

    print *
    print *, "Vacuum polarization estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 4a20)', &
        '<VP>', 'Reference', 'Diff (relative)', 'Tolerance'
    do k = 1, NW
        orb = orbitals(k)
        vp_value = qedvp_kl(k, k)
        print '(i2,":",i2,a2,4es20.10)', &
            k, NP(k), kappa_to_string(NAK(k)), &
            vp_value, orb%vp, reldiff(vp_value, orb%vp), orb%vp_tol
        call test_isequal(success, "VP", vp_value, orb%vp, orb%vp_tol)
    end do

    ! Test against the legacy routines:
    allocate(vp_matrix(NW, NW))
    call qedvp(vp_matrix)
    call qedvp_legacy_init

    do k = 1, NW
        do l = 1, NW
            ! VPINT and VPINTF do not care about the angular momentum symmetry, so we need to
            ! do that check here manually.
            if(NAK(k) /= NAK(l)) then
                label = "NAK(" // itoa(k) // ") /= NAK(" // itoa(l) // ")"
                call test_isequal_atol(success, label, vp_matrix(k, l), 0.0_dp, 0.0_dp)
                cycle
            end if
            label = "VPINT(" // itoa(k) // ", " // itoa(l) // ")"
            call vpint_safe(k, l, vp_value)
            if(vp_matrix(k, l) == 0.0_dp) then
                call test_isequal_atol(success, label, vp_value, vp_matrix(k, l), atol)
            else
                call test_isequal(success, label, vp_value, vp_matrix(k, l), rtol)
            end if
        end do
    end do

    if(.not.success) then
        print *, "hydrogenic_vp: Tests failed."
        stop 1
    end if

contains


    !> Calls the `VPINT` routine, but ensures that the input `k` and `l` arguments do not
    !! get swapped.
    subroutine vpint_safe(k, l, result)
        use vpint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call VPINT(k, l, result)
    end

    subroutine allocate_hydrogenic_orbitals
        use parameter_def
        use memory_man
        use orb_C
        use wave_C

        type(orbital_reference) :: orb
        integer :: k

        NW = size(orbitals)
        print *, "Allocating hydrogenic orbitals: NW=", NW
        call alloc(PF, NNNP, NW, 'PF', 'LODRWF')
        call alloc(QF, NNNP, NW, 'QF', 'LODRWF')

        do k = 1, NW
            orb = orbitals(k)
            call populate_hydrogenic(k, orb%orbital_n, orb%orbital_kappa)
        enddo
    end subroutine allocate_hydrogenic_orbitals

    subroutine populate_hydrogenic(idx, orbital_n, orbital_kappa)
        use def_C
        use wave_C
        use dcbsrw_I

        integer, intent(in) :: idx, orbital_n, orbital_kappa
        real(real64) :: energy, dcwf_RG0
        integer :: i

        NP(idx) = orbital_n
        NAK(idx) = orbital_kappa

        ! Populates the first MF(idx) points in PF(:, idx) and QF(:, idx) with the corresponding
        ! DC wavefunction for n / kappa / Z value.
        call dcbsrw(orbital_n, orbital_kappa, Z, energy, dcwf_RG0, PF(:, idx), QF(:, idx), MF(idx))
    end subroutine populate_hydrogenic

end program hydrogenic_vp
