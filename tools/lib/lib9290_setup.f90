!> Contains routines that can be used to set up the different global state
!! required by the lib9290 routines.
module grasptest_lib9290_setup
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    !> Call all the setup routines to set up parts of the lib9290 global state
    !! for a point nucleus with charge `nuclear_z` and mass `nuclear_mass` (amu).
    !! Nuclear mass is used for mass shifts -- the charge distribution is still
    !! assumed to be a point.
    subroutine setup(nuclear_z, nuclear_mass)
        use grasp_rciqed_lib9290_init

        real(real64), intent(in) :: nuclear_z, nuclear_mass

        call lib9290_init_constants
        call lib9290_init_grid(nuclear_z)
        call lib9290_init_nucleus(nuclear_z)
        call lib9290_init_nucleus_mass(nuclear_mass)
    end subroutine setup

    function kappa_to_string(kappa)
        integer, intent(in) :: kappa
        character, parameter, dimension(*) :: spec_notation = ['s','p','d','f','g','h','i','k','l','m','n']
        character(2) :: kappa_to_string

        if(kappa > 0) then
            kappa_to_string(1:1) = spec_notation(abs(kappa) + 1)
            kappa_to_string(2:2) = '-'
        else
            kappa_to_string(1:1) = spec_notation(abs(kappa))
            kappa_to_string(2:2) = ' '
        end if
    end function kappa_to_string

end module grasptest_lib9290_setup
