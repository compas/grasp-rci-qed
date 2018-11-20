set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_library(rci-qed_matrixelements STATIC
    cimatrixelements.f90
)
setup_fortran_modules(rci-qed_matrixelements)
target_link_libraries_Fortran(rci-qed_matrixelements PRIVATE
    mod 9290 rang90
    grasp rci_rci rci_breit rci_qed rci_mass_shifts
)
