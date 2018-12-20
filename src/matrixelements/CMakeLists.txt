# Warn if implicit interfaces are boung used.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_library(rci-qed_matrixelements STATIC
    cimatrixelements.f90
    rcisettings.f90
    # Legacy code
    iabint.f90 iabint_I.f90
    rkintc.f90 rkintc_I.f90
)
setup_fortran_modules(rci-qed_matrixelements)
target_link_libraries_Fortran(rci-qed_matrixelements PRIVATE
    mod 9290 rang90
    grasp rci_rci rci_breit rci_qed rci_mass_shifts toml
)
