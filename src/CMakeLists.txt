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

# rci-qed.matrixelements
add_executable(rci-qed.pt
    pt.f90
)
target_link_libraries_Fortran(rci-qed.pt PUBLIC
    mod 9290 rang90
    rciqed-grasp rciqed-rci rciqed-breit rciqed-qed rciqed-mass_shifts rciqed-matrixelements
)
install(TARGETS rci-qed.pt DESTINATION bin/)
