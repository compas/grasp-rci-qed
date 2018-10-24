add_subdirectory("breit")
add_subdirectory("qed")
add_subdirectory("mass_shifts")
add_subdirectory("rci")

add_executable(rci-qed
    rci-qed.f90
)
target_link_libraries_Fortran(rci-qed PUBLIC
    mod 9290 mpi90
    rci
)
install(TARGETS rci-qed DESTINATION bin/)
