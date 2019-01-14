!> Functions and routines related to nuclear models etc.
module grasp_rciqed_nucleus
    implicit none

contains

    !> Estimates the RMS of a Fermi nucleus, given the `a` and `c` parameters.
    function fermi_rms(a, c)
        use grasp_rciqed_kinds, only: real64
        use estrms_I, only: ESTRMS

        real(real64), intent(in) :: a, c
        real(real64) :: fermi_rms

        fermi_rms = ESTRMS(a, c)
    end function fermi_rms

end module grasp_rciqed_nucleus
