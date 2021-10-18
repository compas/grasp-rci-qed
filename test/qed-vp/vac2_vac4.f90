program qedvp_vac2_vac4_test
    use grasp_rciqed_lib9290_init
    ! GRASP global variables:
    use parameter_def, only: NNNP
    use grid_C, only: R
    use grasp_rciqed_qed_vp, only: qedvp_init, ZDIST, vp_vac2, vp_vac4
    ! Testing libraries
    use grasptest_testing
    use grasptest_utilities
    implicit none
    ! Settings and reference values that are initialized by one of the init routines below:
    integer, dimension(:), allocatable :: idxs
    real(real64), dimension(:), allocatable :: ref_r, ref_zdist, ref_vac2, ref_vac2p4
    ! Other local variables, including the success state variable:
    character(len=:), allocatable :: testset
    integer :: i
    logical :: success = .true.

    ! Calls SETCON and SETMC:
    call lib9290_init_constants

    ! Initialize setting and reference values, based on the test case:
    if(.not.get_command_argument_allocating(1, testset)) then
        print '("ERROR: unable to extract the command line argument.")'
        error stop
    endif
    if(testset == "pnc_118") then
        call init_refs_pnc_118
    elseif(testset == "fnc_92") then
        call init_refs_fnc_92
    else
        print '("ERROR: invalid testset ",a)', testset
        error stop
    endif

    ! Testing to make sure that the grid is the same. The point is that we are testing
    ! VP values at specific grid points, so
    call test_array_subset("R", idxs, R, ref_r)

    ! Initialize the VP arrays and test their contents:
    call qedvp_init
    call test_array_subset("ZDIST", idxs, ZDIST, ref_zdist)
    call test_array_subset("[VAC2]", idxs, vp_vac2, ref_vac2)
    call test_array_subset("[VAC2+VAC4]TB", idxs, vp_vac2 + vp_vac4, ref_vac2p4)

    if(.not.success) then
        print *, "qedvp_vac2_vac4_test: Tests failed."
        stop 1
    end if

contains

    ! This routine loops over the `TB` array and tests against the values in `refs`. Only
    ! the indices in `idxs` are tested.
    !
    ! Modified the program global variable `success`.
    subroutine test_array_subset(name, idxs, values, refs)
        use tatb_C, only: TB
        use grasptest_testing
        ! Tolerances:
        real(real64), parameter :: atol = 1e-15_dp, rtol = 1e-15_dp
        ! Arguments:
        character(*), intent(in) :: name
        integer, dimension(:), intent(in) :: idxs
        real(real64), dimension(:), intent(in) :: values, refs
        ! Local variables
        character(len=:), allocatable :: loop_label
        integer :: i

        ! Just to make sure that the input arrays are of the same size:
        call test_isequal(success, name//":size(idxs) == size(refs)", size(idxs), size(refs))

        do i = 1, size(idxs)
            if(idxs(i) > size(values)) then
                print '("  Index too larger in ",a,": accessing at ",I0," for size(values)=",I0,")")', &
                    name, idxs(i), size(values)
                success = .false.
                cycle
            endif
            loop_label = name // "(" // itoa(idxs(i)) // ")"
            if(refs(i) == 0.0_dp) then
                call test_isequal_atol(success, loop_label, values(idxs(i)), 0.0_dp, rtol)
            else
                call test_isequal(success, loop_label, values(idxs(i)), refs(i), atol)
            endif
        enddo
    end

    ! Reference data:

    subroutine init_refs_pnc_118
        real(real64) :: z = 118.0_dp
        print *, "Test case: PNC with Z=118"
        idxs = (/ 1, 2, 10, 100, 200, 300, 400 /)
        ref_r = (/ &
            0.0_dp, &
            8.6900163349193419e-10_dp, &
            9.6324099235621998e-9_dp, &
            2.3758468461267459e-6_dp, &
            3.5510546409794880e-4_dp, &
            5.2704822269483954e-2_dp, &
            7.8220916714821236e0_dp &
        /)
        ! For PNC, the ZDIST arrays is just all zeroes
        ref_zdist = (/ ( 0.0_dp, i=1,7 ) /)
        ref_vac2 = (/ &
            0.0_dp, &
            -3055887846.7541680_dp, &
            -230057769.00964698_dp, &
            -509160.07519363466_dp, &
            -886.88203951022888_dp, &
            -5.1178739027808850e-8_dp,  &
            0.0_dp &
        /)
        ref_vac2p4 = (/ &
            0.0_dp, &
            -3131285190.0396390_dp, &
            -234896387.46195826_dp, &
            -515741.12339147180_dp, &
            -893.25091896918639_dp, &
            -5.2166645062196083e-8_dp, &
            0.0_dp &
        /)

        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(z)
        ! Set NPARM = 0 and calls NUCPOT:
        call lib9290_init_nucleus_pnc(z)
    end

    subroutine init_refs_fnc_92
        real(real64) :: z = 92.0_dp
        real(real64) :: a = 9.890591365741376e-6_dp, c = 0.00014178890920967197_dp
        print *, "Test case: FNC with Z=92"
        idxs = (/ 1, 2, 10, 120, 195, 305, 380 /)
        ref_r = (/ &
            0.0_dp, &
            1.1145890516526981e-9_dp, &
            1.2354612728047168e-8_dp, &
            8.3207247621981630e-6_dp, &
            3.5470885213077558e-4_dp, &
            8.6799692069805201e-2_dp, &
            3.6908177267126243_dp &
        /)
        ref_zdist = (/ &
            0.0_dp, &
            1.8727187006528125e-4_dp, &
            3.0967397116672744e-3_dp, &
            510.33579825161564_dp, &
            4.1389190360958001e-4_dp, &
            0.0_dp, &
            0.0_dp &
        /)
        ref_vac2 = (/ &
            0.0_dp, &
            -4438.7607282950403_dp, &
            -4438.7607150539343_dp, &
            -4432.9079728155803_dp, &
            -700.00744353361927_dp, &
            -1.0812577738386876e-12_dp, &
            0.0_dp &
        /)
        ref_vac2p4 = (/ &
            0.0_dp, &
            -4476.2380030732465_dp, &
            -4476.2379897101273_dp, &
            -4470.3352280316485_dp, &
            -705.05629777606441_dp, &
            -1.1060619437733789e-12_dp, &
            0.0_dp &
        /)

        ! Sets up the necessary values and calls SETQIC and RADGRD:
        call lib9290_init_grid(z)
        ! Set NPARM to 2
        call lib9290_init_nucleus_fnc(z, a, c)
    end

end program
