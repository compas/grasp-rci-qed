program qedvp_vac2_vac4_test
    use grasp_rciqed_kinds, only: real64, dp
    use grasp_rciqed_lib9290_init
    ! GRASP global variables:
    use parameter_def, only: NNNP
    use grid_C, only: R, RP
    use tatb_C, only: MTP, TB
    use ncdist_C, only: ZDIST
    ! Routines:
    use ncharg_I
    use vac2_I
    use vac4_I
    ! Testing libraries
    use grasptest_testing
    implicit none
    ! Settings and reference values that are initialized by one of the init routines below:
    real(real64) :: nuclear_Z
    integer, dimension(:), allocatable :: idxs
    real(real64), dimension(:), allocatable :: ref_r, ref_zdist, ref_vac2, ref_vac2p4
    ! Other local variables, including the success state variable:
    integer :: i
    logical :: success = .true.

    ! Initialize setting and reference values, based on the test case:
    call init_refs_pnc_1

    ! Calls SETCON and SETMC:
    call lib9290_init_constants
    ! Sets up the necessary values and calls SETQIC and RADGRD:
    call lib9290_init_grid(nuclear_z)
    ! Set NPARM = 0 and calls NUCPOT:
    call lib9290_init_nucleus_pnc(nuclear_z)

    ! Testing to make sure that the grid is the same. The point is that we are testing
    ! VP values at specific grid points, so
    call test_array_subset("R", idxs, R, ref_r)

    ! Initialize the nuclear charge distribution variable ZDIST. It shoud be all zeroes
    ! if it's PNC:
    call NCHARG
    ! This "re-definition" of ZDIST is in VACPOL.
    ZDIST(:MTP) = ZDIST(:MTP)*R(:MTP)*RP(:MTP)
    call test_array_subset("ZDIST", idxs, ZDIST, ref_zdist)

    ! VAC2, which should be calculating the Uehling potential, populates the TB array:
    call VAC2
    call test_array_subset("[VAC2]TB", idxs, TB, ref_vac2)

    ! VAC4, which should be calculating the Källen-Sabry potential, adds its potential
    ! to the TB array, on top of Uehling. So the values we are testing here are U + KS.
    call VAC4
    call test_array_subset("[VAC2+VAC4]TB", idxs, TB, ref_vac2p4)

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
                continue
            endif
            loop_label = name // "(" // itoa(idxs(i)) // ")"
            if(refs(i) == 0.0_dp) then
                call test_isequal_atol(success, loop_label, values(idxs(i)), 0.0_dp, rtol)
            else
                call test_isequal(success, loop_label, values(idxs(i)), refs(i), atol)
            endif
        enddo
    end subroutine test_array_subset

    ! Borrowed from https://stackoverflow.com/a/31028207
    function itoa(i)
        character(len=:), allocatable :: itoa
        integer,intent(in) :: i
        character(range(i)+2) :: tmp
        write(tmp,'(i0)') i
        itoa = trim(tmp)
    end function itoa

    subroutine init_refs_pnc_1
        print *, "Test case: PNC with Z=118"
        nuclear_Z = 118.0_dp
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
    end subroutine init_refs_pnc_1

end program qedvp_vac2_vac4_test
