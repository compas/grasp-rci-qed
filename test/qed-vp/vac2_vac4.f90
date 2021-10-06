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

    real(real64) :: nuclear_Z = 118.0_dp
    !real(real64) :: vac2_arr(NNNP), vac4_arr(NNNP)
    integer, dimension(7) :: idxs = (/ 1, 2, 10, 100, 200, 300, 400 /)
    real(real64), dimension(7) :: refvals
    integer :: i
    logical :: success = .true.

    ! Calls SETCON and SETMC:
    call lib9290_init_constants
    ! Sets up the necessary values and calls SETQIC and RADGRD:
    call lib9290_init_grid(nuclear_z)
    ! Set NPARM = 0 and calls NUCPOT:
    call lib9290_init_nucleus_pnc(nuclear_z)

    ! Testing to make sure that the grid is the same. The point is that we are testing
    ! VP values at specific grid points, so
    call test_isequal_atol(success, "R(1)", R(1), 0.0000000000000000e0_dp, 1e-15_dp)
    refvals(2:7) = (/ &
        8.6900163349193419e-10_dp, 9.6324099235621998e-9_dp, 2.3758468461267459e-6_dp, &
        3.5510546409794880e-4_dp, 5.2704822269483954e-2_dp, 7.8220916714821236e0_dp &
    /)
    do i = 2, size(idxs)
        call test_isequal_atol(success, "R(i)", R(idxs(i)), refvals(i), 1e-15_dp)
    enddo

    ! Initialize the nuclear charge distribution variable ZDIST. It shoud be all zeroes
    ! if it's PNC:
    call NCHARG
    ! This "re-definition" of ZDIST is in VACPOL.
    ZDIST(:MTP) = ZDIST(:MTP)*R(:MTP)*RP(:MTP)
    ! For PNC, the ZDIST arrays is just all zeroes
    do i = 1, size(idxs)
        call test_isequal_atol(success, "ZDIST(i)", ZDIST(i), 0.0000000000000000e0_dp, 1e-15_dp)
    enddo

    ! VAC2, which should be calculating the Uehling potential, populates the TB array:
    call VAC2
    refvals(:) = (/ &
        0.0_dp, -3055887846.7541680_dp, -230057769.00964698_dp, -509160.07519363466_dp, &
        -886.88203951022888_dp, -5.1178739027808850e-8_dp, 0.0_dp &
    /)
    do i = 1, size(idxs)
        if(refvals(i) == 0.0_dp) then
            call test_isequal_atol(success, "TB(i)", TB(idxs(i)), 0.0_dp, 1e-15_dp)
        else
            call test_isequal(success, "TB(i)", TB(idxs(i)), refvals(i), 1e-15_dp)
        endif
    enddo

    ! VAC4, which should be calculating the Källen-Sabry potential, adds its potential
    ! to the TB array, on top of Uehling. So the values we are testing here are U + KS.
    call VAC4
    call test_isequal_atol(success, "TB(1)", TB(1), 0.0000000000000000e0_dp, 1e-15_dp)
    refvals(:) = (/ &
        0.0_dp, -3131285190.0396390_dp, -234896387.46195826_dp, -515741.12339147180_dp, &
        -893.25091896918639_dp, -5.2166645062196083e-8_dp, 0.0_dp &
    /)
    do i = 1, size(idxs)
        if(refvals(i) == 0.0_dp) then
            call test_isequal_atol(success, "TB(i)", TB(idxs(i)), 0.0_dp, 1e-15_dp)
        else
            call test_isequal(success, "TB(i)", TB(idxs(i)), refvals(i), 1e-15_dp)
        endif
    enddo

    if(.not.success) then
        print *, "qedvp_vac2_vac4_test: Tests failed."
        stop 1
    end if

end program qedvp_vac2_vac4_test
