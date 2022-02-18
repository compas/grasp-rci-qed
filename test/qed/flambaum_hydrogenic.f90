program qed_flambaum_hydrogenic_test
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNN1, NNNP
    use orb_C, only: NW, NP, NAK
    use grasp_rciqed_qed_flambaum
    use grasp_rciqed_lib9290_init
    use grasptest_lib9290_setup, only: kappa_to_string
    use grasptest_testing
    implicit none

    type orbital_reference
        integer :: orbital_n, orbital_kappa
        real(real64) :: se, se_tol ! self-energy and its tolerance
    end type orbital_reference

    ! This reference data was simply calculated using the same code it is testing.
    ! I.e. it is not from a separate reference or an independent calculation and
    ! therefore one should not assume that these numbers are necessarily correct.
    ! They are here to catch any change to the diagonal element values.
    type(orbital_reference), parameter, dimension(*) :: orbitals = (/ &
        orbital_reference(1, -1,  4.4735233e-02, 1e-7), & ! 1s
        orbital_reference(2, -1,  6.0076139e-03, 1e-7), & ! 2s
        orbital_reference(2,  1, -9.4509135e-05, 1e-7), & ! 2p-
        orbital_reference(2, -2,  2.9484379e-04, 1e-7), & ! 2p
        orbital_reference(3,  2, -1.6035006e-05, 1e-7), & ! 3d-
        orbital_reference(3, -3,  2.4125943e-05, 1e-7), & ! 3d
        orbital_reference(4,  3, -4.5917906e-06, 1e-7), & ! 4f-
        orbital_reference(5, -4,  2.0465297e-06, 1e-7) &  ! 5f
    /)
    type(orbital_reference) :: orb

    real(real64), parameter :: nuclear_z = 18.0_dp

    logical :: success = .true.
    integer :: n, k, l
    real(real64) :: qedse, phi_l, phi_f, phi_g
    real(real64) :: qedse_pyykkoe, pyyk, mohr, qedse_mohr_ukw

    call lib9290_init_constants
    call lib9290_init_grid(nuclear_z)
    call init_nucleus ! Set up a Fermi nucleus
    call allocate_orbitals

    PRINT *, size(orbitals)

    print *
    print *, "Flambaum SE estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 4a15)', &
        'Flambaum', 'Phi_l', 'Phi_f', 'Phi_g'
    do k = 1, NW
        orb = orbitals(k)
        qedse = qedse_flambaum(k, k, phi_l, phi_f, phi_g)
        print '(i2,":",i2,a2,4es15.7)', &
            k, NP(k), kappa_to_string(NAK(k)), &
            qedse, phi_l, phi_f, phi_g
        call test_isequal(success, "self-energy", qedse, orb%se, orb%se_tol)
    end do

    if(.not.success) then
        print *, "qed_flambaum_hydrogenic_test: Tests failed."
        stop 1
    end if

contains

    subroutine init_nucleus
        use def_C, only: CVAC, C, PI, TENMAX, EXPMAX, EXPMIN, PRECIS, Z
        use npar_C, only: NPARM, PARM
        use nucpot_I

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        print *, TENMAX,EXPMAX,EXPMIN,PRECIS
        print *, CVAC, PI

        Z = nuclear_z

        NPARM = 2
        PARM(1) = 6.8839456274865651D-005
        PARM(2) = 9.8905913700962641D-006

        call NUCPOT
    end subroutine init_nucleus

    subroutine allocate_orbitals
        use grasptest_lib9290_hydrogenic

        type(orbital_definition), dimension(size(orbitals)) :: orbital_defs

        do k = 1, size(orbitals)
            orbital_defs(k)%n = orbitals(k)%orbital_n
            orbital_defs(k)%kappa = orbitals(k)%orbital_kappa
        enddo

        call allocate_hydrogenic_orbitals(orbital_defs)

    end subroutine allocate_orbitals

end program qed_flambaum_hydrogenic_test
