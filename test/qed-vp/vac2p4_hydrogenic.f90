program qed_vp_hydrogenic_test
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNNW, NNN1, NNNP
    use orb_C, only: NW, NP, NAK
    use grasp_rciqed_qed, only: init_vacuum_polarization
    use grasp_rciqed_qed_vp
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
        orbital_reference(1, -1, -3.1701693522842569e-3_dp, 1e-15), & ! 1s
        orbital_reference(2, -1, -4.0115744750909485e-4_dp, 1e-15), & ! 2s
        orbital_reference(2,  1, -1.8052467524874912e-6_dp, 1e-15), & ! 2p-
        orbital_reference(2, -2, -3.4401849263894386e-7_dp, 1e-15), & ! 2p
        orbital_reference(3,  2, -1.1242875292363751e-10_dp, 1e-15), & ! 3d-
        orbital_reference(3, -3, -3.3481540980321995e-11_dp, 1e-15), & ! 3d
        orbital_reference(4,  3, -6.9334609072624289e-15_dp, 1e-15), & ! 4f-
        orbital_reference(5, -4, -2.0900526959263732e-15_dp, 1e-15)  & ! 5f
    /)
    type(orbital_reference) :: orb

    real(real64), parameter :: nuclear_z = 18.0_dp

    logical :: success = .true.
    integer :: n, k, l
    real(real64), dimension(:, :), allocatable :: qedvpkl
    real(real64) :: vpint

    call lib9290_init_constants
    call lib9290_init_grid(nuclear_z)
    call init_nucleus ! Set up a Fermi nucleus
    call allocate_orbitals

    ! Populate a matrix with the matrix elements of the QED VP:
    call init_vacuum_polarization
    allocate(qedvpkl(NW, NW))
    call qedvp(qedvpkl)

    print *
    print *, "QED VP estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 1a15)', 'QED VP'
    do k = 1, NW
        orb = orbitals(k)
        print '(i2,":",i2,a2,es15.7)', &
            k, NP(k), kappa_to_string(NAK(k)), qedvpkl(k, k)
        call test_isequal(success, "qedvp", qedvpkl(k, k), orb%se, orb%se_tol)
    end do

    ! Here we test against the old GRASP routines, specifically the values given by the
    ! VPINT routine. However, one thing to note is that the old routines ignore the
    ! angular momentum quantum number. This is likely because the angular coefficient
    ! calculation already takes that into account. But we have to implement the check
    ! ourselves here.
    print *
    print *, "QED VP comparisons with old VP routines:"
    print '(11(" "), 2a15)', 'VPINT', 'qedvp'
    do k = 1, NW
        do l = 1, NW
            if(NAK(k) /= NAK(l)) then
                cycle
            endif
            call vpint_safe(k, l, vpint)
            print '(i2,a2,":",i2,a2,": ",2es15.7)', &
                NP(k), kappa_to_string(NAK(k)), NP(k), kappa_to_string(NAK(k)), &
                vpint, qedvpkl(k, l)
            if(qedvpkl(k, l) == 0.0_dp) then
                call test_isequal_atol(success, "qedvp~vpint", vpint, qedvpkl(k, l), 1e-15_dp)
            else
                call test_isequal(success, "qedvp~vpint", vpint, qedvpkl(k, l), 1e-15_dp)
            endif
        end do
    end do

    if(.not.success) then
        print *, "qed_vp_hydrogenic_test: Tests failed."
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

    !> Calls the `VPINT` routine, but ensures that the input `k` and `l`
    !! arguments do not get swapped.
    subroutine vpint_safe(k, l, result)
        use grasp_rciqed_kinds, only: real64
        use vpint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call VPINT(k, l, result)
    end

end program
