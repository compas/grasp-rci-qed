# Set up a variable that points to the root of the grasp-rci-qed repository
get_filename_component(GRASPRCI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
message("* GRASPRCI_ROOT: ${GRASPRCI_ROOT}")

add_subdirectory("libgrasp")
add_subdirectory("breit")
add_subdirectory("qed")
add_subdirectory("mass_shifts")
add_subdirectory("rci")
add_subdirectory("matrixelements")
add_subdirectory("toml")

# rci-qed is the main program in this module
add_executable(rci-qed
    rci-qed.f90
)
target_link_libraries_Fortran(rci-qed PUBLIC
    mod 9290 mpi90
    rciqed-rci rciqed-matrixelements
)
install(TARGETS rci-qed DESTINATION bin/)

# rci-qed.pt
add_executable(rci-qed.pt
    pt.f90
)
target_link_libraries_Fortran(rci-qed.pt PUBLIC
    mod 9290 rang90
    rciqed-grasp rciqed-rci rciqed-breit rciqed-qed rciqed-mass_shifts rciqed-matrixelements
)
install(TARGETS rci-qed.pt DESTINATION bin/)

# rci-qed.orbitals
add_executable(rci-qed.orbitals
    orbitals.F90
)
target_link_libraries_Fortran(rci-qed.orbitals PUBLIC
    mod 9290 rang90
    rciqed-grasp rciqed-rci rciqed-breit rciqed-qed rciqed-mass_shifts rciqed-matrixelements
)
install(TARGETS rci-qed.orbitals DESTINATION bin/)

add_subdirectory("libgrasp-rci")
