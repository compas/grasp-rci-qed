!> Creates a `.w` file with hydrogenic wavefunctions.
program exporthydrogenic
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

end program exporthydrogenic
