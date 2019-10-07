set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface")

add_library(grasp-rci SHARED
    libgrasp-rci.f90
    libgrasprci.f90 grasp_c_interface.f90
)
target_link_libraries_Fortran(grasp-rci PUBLIC
    mod 9290 rang90
    rciqed-grasp rciqed-rci rciqed-breit rciqed-qed rciqed-mass_shifts rciqed-matrixelements
)
install(TARGETS grasp-rci DESTINATION lib/)
install(TARGETS grasp-rci DESTINATION "${GRASPRCI_ROOT}/lib/")
