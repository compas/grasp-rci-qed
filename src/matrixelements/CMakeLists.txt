set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_library(rci-qed_matrixelements STATIC
    cimatrixelements.f90
    g2k_lib92.f90
    g2k_librci.f90
    g2k_csls.F90
    eval_matrixelement.f90
    # Legacy code
    getmixblock.f90 getmixblock_I.f90
)
setup_fortran_modules(rci-qed_matrixelements)
target_link_libraries_Fortran(rci-qed_matrixelements PRIVATE
    mod 9290 rang90
    grasp rci_rci rci_breit rci_qed rci_mass_shifts
)
