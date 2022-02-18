program qed_vp_hydrogenic_test
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use parameter_def, only: NNNW, NNN1, NNNP
    use orb_C, only: NW, NP, NAK
    use grasp_rciqed_qed_vp, only: qedvp_init
    use grasp_rciqed_qed_vp
    use grasp_rciqed_lib9290_init, only: lib9290_init_constants
    use grasptest_lib9290_setup, only: kappa_to_string
    use grasptest_testing
    use grasptest_utilities
    use grasptest_qedvp_legacy, only: qedvp_legacy_init
    implicit none

    type orbital_reference
        integer :: orbital_n, orbital_kappa
        real(real64) :: x, rtol ! reference value and its tolerance
    end type

    ! Settings and reference values that are initialized by one of the init routines below:
    type(orbital_reference), dimension(:), allocatable :: orbitals
    ! Other local variables, including the success state variable:
    type(orbital_reference) :: orb
    character(len=:), allocatable :: testset
    integer :: n, k, l
    real(real64), dimension(:, :), allocatable :: qedvpkl
    real(real64) :: vpint
    logical :: success = .true.

    ! Calls SETCON and SETMC:
    call lib9290_init_constants

    if(.not.get_command_argument_allocating(1, testset)) then
        print '("ERROR: unable to extract the command line argument.")'
        error stop
    endif

    call lib9290_init_constants
    if(testset == "pnc_92") then
        call init_reference_pnc92
    elseif(testset == "fnc_18") then
        call init_reference_fnc18
    else
        print '("ERROR: invalid testset ",a)', testset
        error stop
    endif

    ! Allocating the hydrogenic orbitals:
    call allocate_orbitals

    ! Populate a matrix with the matrix elements of the QED VP:
    call qedvp_init
    allocate(qedvpkl(NW, NW))
    call qedvp(qedvpkl)

    print *
    print *, "QED VP estimates for hydrogenic wavefunctions (DCBSRW routine)"
    print '(7(" "), 1a15)', 'QED VP'
    do k = 1, NW
        orb = orbitals(k)
        print '(i2,":",i2,a2,es23.15)', &
            k, NP(k), kappa_to_string(NAK(k)), qedvpkl(k, k)
        call test_isequal(success, "qedvp", qedvpkl(k, k), orb%x, orb%rtol)
    end do

    ! Here we test against the old GRASP routines, specifically the values given by the
    ! VPINT routine. However, one thing to note is that the old routines ignore the
    ! angular momentum quantum number. This is likely because the angular coefficient
    ! calculation already takes that into account. But we have to implement the check
    ! ourselves here.
    print *
    print *, "QED VP comparisons with old VP routines:"
    call qedvp_legacy_init
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
        use vpint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call VPINT(k, l, result)
    end

    ! Reference data:
    !
    ! This reference data was simply calculated using the same code it is testing.
    ! I.e. it is not from a separate reference or an independent calculation and
    ! therefore one should not assume that these numbers are necessarily correct.
    ! They are here to catch any change to the diagonal element values.
    subroutine init_reference_fnc18
        use grasp_rciqed_lib9290_init
        real(real64) :: z = 18.0_dp
        real(real64) :: a = 9.8905913700962641e-6_dp, c = 6.8839456274865651e-5_dp
        print *, "Test case: FNC with Z=18"
        orbitals = (/ &
            orbital_reference(1, -1, -3.1701693522842569e-3_dp, 1e-15), & ! 1s
            orbital_reference(2, -1, -4.0115744750909485e-4_dp, 1e-15), & ! 2s
            orbital_reference(2,  1, -1.8052467524874912e-6_dp, 1e-15), & ! 2p-
            orbital_reference(2, -2, -3.4401849263894386e-7_dp, 1e-15), & ! 2p
            orbital_reference(3,  2, -1.1242875292363751e-10_dp, 1e-15), & ! 3d-
            orbital_reference(3, -3, -3.3481540980321995e-11_dp, 1e-15), & ! 3d
            orbital_reference(4,  3, -6.9334609072624289e-15_dp, 1e-15), & ! 4f-
            orbital_reference(5, -4, -2.0900526959263732e-15_dp, 1e-15)  & ! 5f
        /)

        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(z)
        ! Set NPARM to 2
        call lib9290_init_nucleus_fnc(z, a, c)
    end

    subroutine init_reference_pnc92
        use grasp_rciqed_lib9290_init
        real(real64) :: z = 92.0_dp
        print *, "Test case: PNC with Z=92"
        orbitals = (/ &
            orbital_reference(1, -1, -3.628948611886053e0_dp, 1e-15), & ! 1s
            orbital_reference(2, -1, -6.409362427230237e-1_dp, 1e-15), & ! 2s
            orbital_reference(2,  1, -1.109676026997344e-1_dp, 1e-15), & ! 2p-
            orbital_reference(2, -2, -4.699013803610754e-3_dp, 1e-15), & ! 2p
            orbital_reference(3,  2, -5.336765778958868e-5_dp, 1e-15), & ! 3d-
            orbital_reference(3, -3, -9.583019021422298e-6_dp, 1e-15), & ! 3d
            orbital_reference(4,  3, -6.208836467501664e-8_dp, 1e-15), & ! 4f-
            orbital_reference(5, -4, -1.416895693533666e-8_dp, 1e-15)  & ! 5f
        /)

        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(z)
        call lib9290_init_nucleus_pnc(z)
    end

end program
