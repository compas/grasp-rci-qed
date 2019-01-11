!> Contains the routines necessary to populate the orb_C and wave_C modules with
!! hydrogenic orbitals.
!!
!! Uses the DCBSRW routine from lib9290 to numerically solve for the wavefunction.
module grasptest_lib9290_hydrogenic
    implicit none

    type orbital_definition
        integer :: n, kappa
    end type orbital_definition

contains

    subroutine allocate_hydrogenic_orbitals(orbitals)
        use parameter_def
        use memory_man
        use orb_C
        use wave_C

        type(orbital_definition), dimension(:), intent(in) :: orbitals
        integer :: k

        NW = size(orbitals)
        print *, "Allocating hydrogenic orbitals: NW=", NW
        call alloc(PF, NNNP, NW, 'PF', 'LODRWF')
        call alloc(QF, NNNP, NW, 'QF', 'LODRWF')

        do k = 1, NW
            call populate_hydrogenic(k)
        enddo

    contains

        subroutine populate_hydrogenic(idx)
            use grasp_rciqed_kinds, only: real64
            use def_C
            use wave_C
            use dcbsrw_I

            integer, intent(in) :: idx
            real(real64) :: energy, dcwf_RG0

            NP(idx) = orbitals(k)%n
            NAK(idx) = orbitals(k)%kappa

            ! Populates the first MF(idx) points in PF(:, idx) and QF(:, idx) with the corresponding
            ! DC wavefunction for n / kappa / Z value.
            call dcbsrw(orbitals(k)%n, orbitals(k)%kappa, Z, energy, dcwf_RG0, PF(:, idx), QF(:, idx), MF(idx))

            print *, orbitals(k)%n, orbitals(k)%kappa, Z, energy, dcwf_RG0
        end subroutine populate_hydrogenic

    end subroutine allocate_hydrogenic_orbitals

end module grasptest_lib9290_hydrogenic
