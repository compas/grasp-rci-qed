set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")
add_executable(rci.matrixelements
    matrixelements.f90
    g2k_lib92.f90
    g2k_librci.f90
    g2k_csls.F90
    eval_matrixelement.f90
    # Legacy code
    getmixblock.f90 getmixblock_I.f90
)
setup_fortran_modules(grasp)
target_link_libraries_Fortran(rci.matrixelements PUBLIC
    mod 9290 rang90
    grasp rci_rci rci_breit rci_qed rci_mass_shifts
)
install(TARGETS rci.matrixelements DESTINATION bin/)

# add_executable(rci_analysis
#     rci_analysis.f90
#     g2k_lib92.f90
#     g2k_librci.f90
#     g2k_csls.f90
#     eval_matrixelement.f90
#     # Legacy code
#     #getmixblock.f
# )
# target_link_libraries(rci_analysis
#     92 rang grasp rcicommon
# )
# install(TARGETS rci_analysis DESTINATION bin/)
