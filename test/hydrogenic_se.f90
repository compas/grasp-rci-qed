program hydrogenic_se
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use grasptest_testing, only: test_isequal, reldiff
    use grasptest_lib9290_setup
    use parameter_def, only: NNNW
    use orb_C
    use qed_slfen_I
    implicit none

    type orbital_reference
        integer :: orbital_n, orbital_kappa
        real(real64) :: value, rtol ! self-energy and its relative tolerance
    end type orbital_reference

    ! This reference data was simply calculated using the same code it is testing.
    ! I.e. it is not from a separate reference or an independent calculation and
    ! therefore one should not assume that these numbers are necessarily correct.
    ! They are here to catch any change to the diagonal element values.
    ! The reference values are accurate up to 10^{-8}.
    type(orbital_reference), parameter, dimension(*) :: orbitals = (/ &
        orbital_reference(1, -1,  4.4601329718e-02, 1d-7), & ! 1s
        orbital_reference(2, -1,  5.9901019556e-03, 1d-7), & ! 2s
        orbital_reference(2,  1, -1.5770357162e-04, 1d-7), & ! 2p-
        orbital_reference(2, -2,  2.2820996447e-04, 1d-7), & ! 2p
        orbital_reference(3,  2, -2.0691015722e-05, 1d-7), & ! 3d-
        orbital_reference(3, -3,  1.9385607876e-05, 1d-7), & ! 3d
        orbital_reference(4,  3, -4.3616687658e-06, 1d-7), & ! 4f-
        orbital_reference(5, -4,  2.1599904445e-06, 1d-7)  & ! 5f
    /)
    type(orbital_reference) :: orb

    logical :: success = .true.
    integer :: n, k, l, tmpk, tmpl
    real(real64) :: vp_value

    real(real64), dimension(NNNW) :: slfint

    call setup(18.0_dp, 39.9623831225_dp) ! Ar-40
    call allocate_hydrogenic_orbitals

    call setup_mohr_se

    ! Call out to the library to calculate the SE contributions for each orbital.
    call QED_SLFEN(slfint)

    print *
    print *, "Mohr self-energy estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 2a25, a16, a14)', &
        '<SE>', 'Reference', 'Diff (relative)', 'Tolerance'
    do k = 1, NW
        orb = orbitals(k)
        print '(i2,":",i2,a2,2es25.15,es16.5,es14.5)', &
            k, NP(k), kappa_to_string(NAK(k)), &
            slfint(k), orb%value, reldiff(slfint(k), orb%value), orb%rtol
        call test_isequal(success, "SE", slfint(k), orb%value, orb%rtol)
    end do

    if(.not.success) then
        print *, "hydrogenic_se: Tests failed."
        stop 1
    end if

contains

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

    subroutine setup_mohr_se
        use decide_C, only: LSE
        use qedcut_C, only: NQEDMAX
        LSE = .true.
        NQEDMAX = 8
    end subroutine setup_mohr_se

end program hydrogenic_se
