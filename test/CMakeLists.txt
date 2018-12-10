set(BUILD_BINDIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
# Make sure the test libraries and executables do not end up in build/{lib,bin}/
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test binaries and libraries, and
# abort the compilation if any warning is encountered.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-procedure -Werror=implicit-interface")

# A library of test helper routines etc.
add_library(rci-qed.test STATIC
    lib/testing.f90
)
setup_fortran_modules(rci-qed.test)
target_link_libraries_Fortran(rci-qed.test PUBLIC
    grasp
)

# Script-based tests
add_test(NAME test-integration-rci-nitrogen
    COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/grasp-integration/rci.sh
    "nitrogen"
)
set_property(TEST test-integration-rci-nitrogen PROPERTY ENVIRONMENT
    "GRASP=${GRASP}"
    "GRASP_BUILD_BINDIR=${BUILD_BINDIR}"
)

# Check lib92 routines -- this currently acts as a simple integration test,
# making sure that everything actually got compiled properly. The test checks
# that the QUAD routine from lib92 actually produces reasonable numbers.
add_executable(test.lib92_quad
    lib92_quad.f90
)
target_link_libraries_Fortran(test.lib92_quad PUBLIC mod 9290)
add_test(lib92_quad test.lib92_quad)

# Unit tests for the 1-particle QED operators
add_executable(test.hydrogenic_vp
    hydrogenic_vp.f90
)
target_link_libraries_Fortran(test.hydrogenic_vp PUBLIC
    mod 9290 grasp
    rci_rci rci_qed
    rci-qed.tools rci-qed.test
)
add_test(hydrogenic_vp test.hydrogenic_vp)

add_executable(test.hydrogenic_se
    hydrogenic_se.f90
)
target_link_libraries_Fortran(test.hydrogenic_se PUBLIC
    mod 9290 grasp
    rci_rci rci_qed
    rci-qed.tools rci-qed.test
)
add_test(hydrogenic_se test.hydrogenic_se)


# Unit tests for matrix element calculations
add_executable(test.matrixelements_hydrogenic
    matrixelements_hydrogenic.f90
    lib/testing.f90
)
target_link_libraries_Fortran(test.matrixelements_hydrogenic PUBLIC
    mod 9290 grasp
    rci-qed_matrixelements
    rci_rci rci_breit rci_qed rci_mass_shifts
    rci-qed.tools rci-qed.test
)
add_test(matrixelements_hydrogenic_generate ../tools/rci-qed.exporthydrogenic)
add_test(matrixelements_hydrogenic test.matrixelements_hydrogenic)
set_property(TEST matrixelements_hydrogenic PROPERTY ENVIRONMENT
    "RCIQED_TESTDATA=${CMAKE_CURRENT_SOURCE_DIR}/data"
)
set_property(TEST matrixelements_hydrogenic PROPERTY DEPENDS
    matrixelements_hydrogenic_generate
)

# Integration tests
add_subdirectory("integration/")

# Unit tests for the TOML library
add_executable(test.toml
    toml/tests.f90
)
target_link_libraries_Fortran(test.toml PUBLIC
    toml rci-qed.test grasp
)
add_test(
    NAME toml
    COMMAND test.toml
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/toml"
)
