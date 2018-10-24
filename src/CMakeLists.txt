add_library(grasp STATIC
    libgrasp/kinds.f90
)
setup_fortran_modules(grasp)

add_subdirectory("breit")
add_subdirectory("qed")
add_subdirectory("mass_shifts")
add_subdirectory("rci")

# rci-qed is the main program in this module
add_executable(rci-qed
    rci-qed.f90
)
target_link_libraries_Fortran(rci-qed PUBLIC
    mod 9290 mpi90
    rci
)
install(TARGETS rci-qed DESTINATION bin/)

# Utility programs
add_executable(rci.exportvp
    exportvp.f90
)
target_link_libraries_Fortran(rci.exportvp PUBLIC
    mod 9290 mpi90
    rci rci_qed
)
install(TARGETS rci.exportvp DESTINATION bin/)
