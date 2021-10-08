!> Functions and routines related to nuclear models etc.
module grasp_rciqed_nucleus
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    !> Estimates the RMS of a Fermi nucleus, given the `a` and `c` parameters.
    function fermi_rms(a, c)
        use estrms_I, only: ESTRMS

        real(real64), intent(in) :: a, c
        real(real64) :: fermi_rms

        fermi_rms = ESTRMS(a, c)
    end function fermi_rms

end module grasp_rciqed_nucleus
