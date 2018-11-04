program matrixelements_hydrogenic
    use grasp_kinds, only: real64, dp
    use grasptest_lib9290_setup
    use grasptest_lib9290_hydrogenic
    use orbout_I
    implicit none

    type(orbital_definition), parameter, dimension(*) :: orbitals = (/ &
        orbital_definition(1, -1), & ! 1s
        orbital_definition(2, -1), & ! 2s
        orbital_definition(2,  1), & ! 2p-
        orbital_definition(2, -2), & ! 2p
        orbital_definition(3,  2), & ! 3d-
        orbital_definition(3, -3), & ! 3d
        orbital_definition(4,  3), & ! 4f-
        orbital_definition(5, -4)  & ! 5f
    /)

    logical :: tests_passed = .true.
    integer :: n, k, l, tmpk, tmpl
    real(real64) :: vp_value

    call setup(74.0_dp)
    call allocate_hydrogenic_orbitals(orbitals)
    call orbout("hydrogenic.w")

contains

    function reldiff(a, b)
        use grasp_kinds, only: real64
        real(real64), intent(in) :: a, b
        real(real64) :: reldiff
        reldiff = abs(a-b) / max(abs(a), abs(b))
    end function reldiff

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

end program matrixelements_hydrogenic
