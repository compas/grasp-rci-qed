! module qed_flambaum_hydrogenic_test_common
!     use g2k_parameters
!     implicit real*8 (a-h, o-z)
!
!     private
!
!     pointer (PNTRPF,PF(NNNP,NNNW))
!     pointer (PNTRQF,QF(NNNP,NNNW))
!
!     COMMON/ORB2/NCF,NW,PNTRIQ &
!           /ORB4/NP(NNNW),NAK(NNNW) &
!           /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N &
!           /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
!
!     public :: NW, NP, NAK
!     public :: PNTRPF, PF, PNTRQF, QF, MF
!     public :: N
! end module qed_flambaum_hydrogenic_test_common

program hydrogenic_vp
    !use g2k_parameters, only: real64, dp, NNN1, NNNP
    !use g2ktest_lib92_common, only: setup_constants, setup_grid, setup_nucleus, kappa_to_string
    !use qed_flambaum_hydrogenic_test_common, only: NW, NP, NAK
    use grasp_kinds, only: real64
    use grasptest_lib9290_setup
    use orb_C
    use vpint_I
    implicit none

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
        orbital_reference(1, -1, -3.1701678053e-03, 1d-7), & ! 1s
        orbital_reference(2, -1, -4.0115725286e-04, 1d-7), & ! 2s
        orbital_reference(2,  1, -1.8052454393e-06, 1d-7), & ! 2p-
        orbital_reference(2, -2, -3.4401802024e-07, 1d-7), & ! 2p
        orbital_reference(3,  2, -1.1242863616e-10, 1d-7), & ! 3d-
        orbital_reference(3, -3, -3.3481529462e-11, 1d-7), & ! 3d
        orbital_reference(4,  3, -6.9334598114e-15, 1d-7), & ! 4f-
        orbital_reference(5, -4, -2.0900530055e-15, 1d-7)  & ! 5f
    /)
    type(orbital_reference) :: orb

    logical :: tests_passed = .true.
    integer :: n, k, l, tmpk, tmpl
    real(real64) :: vp_value

    call setup_constants
    call setup_grid
    call setup_nucleus
    call allocate_hydrogenic_orbitals

    call setup_vacuum_polarization

    print *
    print *, "Vacuum polarization estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 4a20)', &
        '<VP>', 'Reference', 'Diff (relative)', 'Tolerance'
    do k = 1, NW
        orb = orbitals(k)
        tmpk = k ! necessary back, because VPINT arguments are set to be INOUT
        call VPINT(tmpk, tmpk, vp_value)
        print '(i2,":",i2,a2,4es20.10)', &
            k, NP(k), kappa_to_string(NAK(k)), &
            vp_value, orb%vp, reldiff(vp_value, orb%vp), orb%vp_tol
        call check_tolerance("VP", vp_value, orb%vp, orb%vp_tol)
    end do

    if(.not.tests_passed) then
        print *, "qed_flambaum_hydrogenic_test: Tests failed."
        stop 1
    end if

contains

    function reldiff(a, b)
        use grasp_kinds, only: real64
        real(real64), intent(in) :: a, b
        real(real64) :: reldiff
        reldiff = abs(a-b) / max(abs(a), abs(b))
    end function reldiff

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
        !use qed_flambaum_hydrogenic_test_common
        !use qedse_flambaum_common, only: Z

        integer, intent(in) :: idx, orbital_n, orbital_kappa
        real(real64) :: energy, dcwf_RG0
        integer :: i

        NP(idx) = orbital_n
        NAK(idx) = orbital_kappa

        ! Populates the first MF(idx) points in PF(:, idx) and QF(:, idx) with the corresponding
        ! DC wavefunction for n / kappa / Z value.
        call dcbsrw(orbital_n, orbital_kappa, Z, energy, dcwf_RG0, PF(:, idx), QF(:, idx), MF(idx))
    end subroutine populate_hydrogenic

    subroutine setup_vacuum_polarization
        use decide_C, only: LVP
        use vpilst_C, only: FRSTVP, NVPI
        use ncdist_C, only: ZDIST
        use tatb_C, only: TB
        use grid_C, only: N, RP
        use ncharg_I
        use vacpol_I

        LVP = .true.
        call ncharg
        call vacpol
        ZDIST(2:N) = TB(2:N)*RP(2:N)
        FRSTVP = .TRUE.
        NVPI = 0
    end subroutine setup_vacuum_polarization

    subroutine check_tolerance(which, a, b, relative_tolerance)
        use grasptest_testing, only: within_tolerance

        character(*), intent(in) :: which
        real(real64), intent(in) :: a, b, relative_tolerance
        real(real64) :: relative_difference

        if (.not.within_tolerance(a, b, relative_tolerance)) then
            relative_difference = abs(a-b) / max(abs(a), abs(b))
            print '("  Test failed: ",a," not within tolerance. Rel.diff: ",es12.5,", tol: ",es12.5)', &
                which, relative_difference, relative_tolerance
            tests_passed = .false.
        endif
    end subroutine check_tolerance

end program hydrogenic_vp
