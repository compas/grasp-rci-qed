program qedvp_funk_funl
    use funk_I
    use funl_I
    use grasptest_testing
    implicit none

    logical :: success = .true.

    ! We test the FUNK and FUNL VP routines, which implement certain expressions from the
    ! reference paper. The reference values were calculated with the routines themselves,
    ! so they have not been verified directly. These tests currently just make sure that
    ! FUNK and FUNL routines do not break when mucking around with the code.

    call test_isequal(success, "FUNK(0.0, 0)", FUNK(0.0_dp, 0), 0.88357293375000001_dp,   1e-15_dp)
    call test_isequal(success, "FUNK(1.0, 0)", FUNK(1.0_dp, 0), 0.11253955800219992_dp,   1e-15_dp)
    call test_isequal(success, "FUNK(2.0, 0)", FUNK(2.0_dp, 0), 2.6012736423913444e-2_dp, 1e-15_dp)
    call test_isequal(success, "FUNK(1.0, 1)", FUNK(1.0_dp, 1), 0.17793357898629991_dp,   1e-15_dp)
    call test_isequal(success, "FUNK(2.0, 1)", FUNK(2.0_dp, 1), 3.5932722757371340e-2_dp, 1e-15_dp)
    call test_isequal(success, "FUNK(1.0, 3)", FUNK(1.0_dp, 3), 0.67560907338809995_dp,   1e-15_dp)
    call test_isequal(success, "FUNK(2.0, 3)", FUNK(2.0_dp, 3), 8.4134587179192935e-2_dp, 1e-15_dp)
    call test_isequal(success, "FUNK(1.0, 5)", FUNK(1.0_dp, 5), 5.7773775413882991_dp,    1e-15_dp)
    call test_isequal(success, "FUNK(2.0, 5)", FUNK(2.0_dp, 5), 0.29521922706427051_dp,   1e-15_dp)

    call test_isequal(success, "FUNL(0.0, 0)", FUNL(0.0_dp, 0), 2.0081880000000001_dp,    1e-15_dp)
    call test_isequal(success, "FUNL(1.0, 0)", FUNL(1.0_dp, 0), 0.31667000000000006_dp,   1e-15_dp)
    call test_isequal(success, "FUNL(2.0, 0)", FUNL(2.0_dp, 0), 9.0680835698371709e-2_dp, 1e-15_dp)
    call test_isequal(success, "FUNL(1.0, 1)", FUNL(1.0_dp, 1), 0.42521900000000007_dp,   1e-15_dp)
    call test_isequal(success, "FUNL(2.0, 1)", FUNL(2.0_dp, 1), 0.10657934231593091_dp,   1e-15_dp)

    if(.not.success) then
        print *, "qedvp_funk_funl: Tests failed."
        stop 1
    end if

end program qedvp_funk_funl
