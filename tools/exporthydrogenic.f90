!> Creates a `.w` file with hydrogenic wavefunctions.
program exporthydrogenic
    use grasp_rciqed_kinds, only: real64, dp
    use grasptest_lib9290_setup
    use grasptest_lib9290_hydrogenic
    use orbout_I
    implicit none

    type(orbital_definition), parameter, dimension(*) :: orbitals = (/ &
        orbital_definition(1, -1), & ! 1s

        orbital_definition(2, -1), & ! 2s
        orbital_definition(2,  1), & ! 2p-
        orbital_definition(2, -2), & ! 2p

        orbital_definition(3, -1), & ! 3s
        orbital_definition(3,  1), & ! 3p-
        orbital_definition(3, -2), & ! 3p
        orbital_definition(3,  2), & ! 3d-
        orbital_definition(3, -3), & ! 3d

        orbital_definition(4, -1), & ! 4s
        orbital_definition(4,  1), & ! 4p-
        orbital_definition(4, -2), & ! 4p
        orbital_definition(4,  2), & ! 4d-
        orbital_definition(4, -3), & ! 4d
        orbital_definition(4,  3), & ! 4f-
        orbital_definition(4, -4), & ! 4f

        orbital_definition(5, -1), & ! 5s
        orbital_definition(5,  1), & ! 5p-
        orbital_definition(5, -2), & ! 5p
        orbital_definition(5,  2), & ! 5d-
        orbital_definition(5, -3), & ! 5d

        orbital_definition(6, -1)  & ! 6s
    /)

    logical :: tests_passed = .true.
    integer :: n, k, l, tmpk, tmpl
    real(real64) :: vp_value

    call setup(74.0_dp, 183.91033628717801_dp)
    call allocate_hydrogenic_orbitals(orbitals)
    call orbout("hydrogenic.w")

end program exporthydrogenic
