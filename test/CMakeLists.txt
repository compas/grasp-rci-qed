set(BUILD_BINDIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
# Make sure the test libraries and executables do not end up in build/{lib,bin}/
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test binaries and libraries, and
# abort the compilation if any warning is encountered.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-procedure -Werror=implicit-interface")

# A library of test helper routines etc.
add_library(rciqed-test STATIC
    lib/utilities.f90
    lib/testing.f90
)
setup_fortran_modules(rciqed-test)
target_link_libraries_Fortran(rciqed-test PUBLIC
    rciqed-grasp
)

# A library of legacy GRASP92 routines, used in testing against old implementations.
add_library(rciqed-test-legacy STATIC
    qed-vp/legacy/grasptest_qedvp_legacy.f90
    # libmod modules:
    qed-vp/legacy/vpilst_C.f90
    # routines:
    qed-vp/legacy/vpint_I.f90 qed-vp/legacy/vpint.f90
    qed-vp/legacy/vpintf_I.f90 qed-vp/legacy/vpintf.f90
)
setup_fortran_modules(rciqed-test-legacy)
target_link_libraries_Fortran(rciqed-test-legacy PUBLIC
    mod 9290 rciqed-qed
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
    mod 9290
    rciqed-grasp rciqed-rci rciqed-qed
    rciqed-tools rciqed-test rciqed-test-legacy
)
add_test(hydrogenic_vp test.hydrogenic_vp)

add_executable(test.hydrogenic_se
    hydrogenic_se.f90
)
target_link_libraries_Fortran(test.hydrogenic_se PUBLIC
    mod 9290
    rciqed-grasp rciqed-rci rciqed-qed
    rciqed-tools rciqed-test
)
add_test(hydrogenic_se test.hydrogenic_se)

# Flambaum QED w/ hydrogenic orbitals
add_executable(test.qed_flambaum_hydrogenic
    qed/flambaum_hydrogenic.f90
)
target_link_libraries_Fortran(test.qed_flambaum_hydrogenic PUBLIC
    mod 9290
    rciqed-grasp rciqed-rci rciqed-qed
    rciqed-tools rciqed-test
)
add_test(qed_flambaum_hydrogenic test.qed_flambaum_hydrogenic)

# QEDMOD QED w/ hydrogenic orbitals
add_executable(test.qed_qedmod_hydrogenic
    qed/qedmod_hydrogenic.f90
)
target_link_libraries_Fortran(test.qed_qedmod_hydrogenic PUBLIC
    mod 9290
    rciqed-grasp rciqed-rci rciqed-qed
    rciqed-tools rciqed-test
)
add_test(qed_qedmod_hydrogenic test.qed_qedmod_hydrogenic)

# Tests for vacuum polarization routines
add_executable(test.qedvp_funk_funl
    qed-vp/funk_funl.f90
)
target_link_libraries_Fortran(test.qedvp_funk_funl PUBLIC
    rciqed-qed rciqed-test
)
add_test(qedvp_funk_funl test.qedvp_funk_funl)

add_executable(test.qedvp_vac2_vac4
    qed-vp/vac2_vac4.f90
)
target_link_libraries_Fortran(test.qedvp_vac2_vac4 PUBLIC
    mod
    rciqed-qed rciqed-test
)
add_test(qedvp_vac2_vac4_pnc_118 test.qedvp_vac2_vac4 "pnc_118")
add_test(qedvp_vac2_vac4_fnc_92 test.qedvp_vac2_vac4 "fnc_92")

add_executable(test.qedvp_hydrogenic
    qed-vp/vac2p4_hydrogenic.f90
)
target_link_libraries_Fortran(test.qedvp_hydrogenic PUBLIC
    mod 9290
    rciqed-grasp rciqed-rci rciqed-qed
    rciqed-tools rciqed-test rciqed-test-legacy
)
add_test(qedvp_hydrogenic_fnc_18 test.qedvp_hydrogenic "fnc_18")
add_test(qedvp_hydrogenic_pnc_92 test.qedvp_hydrogenic "pnc_92")

# Unit tests for matrix element calculations
add_subdirectory("matrixelements/")

# Integration tests
add_subdirectory("integration/")

# Unit tests for the TOML library
add_executable(test.toml
    toml/tests.f90
)
target_link_libraries_Fortran(test.toml PUBLIC
    rciqed-grasp rciqed-toml
    rciqed-test
)
add_test(
    NAME toml
    COMMAND test.toml
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/toml"
)
