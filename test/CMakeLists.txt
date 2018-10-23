# Make sure the test executables end up on build/test/, not build/bin/
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test prorams and libraries
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_test(NAME test-integration-rci-nitrogen
    COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/integration/rci.sh
    "nitrogen"
)
set_property(TEST test-integration-rci-nitrogen PROPERTY ENVIRONMENT
    "GRASP=${GRASP}"
    "RCIQED=$<TARGET_FILE:rci-qed>"
)

# Check lib92 routines -- this currently acts as a simple integration test,
# making sure that everything actually got compiled properly. The test checks
# that the QUAD routine from lib92 actually produces reasonable numbers.
add_executable(test.lib92_quad
    lib92_quad.f90
)
target_link_libraries_Fortran(test.lib92_quad PUBLIC mod 9290)
add_test(lib92_quad test.lib92_quad)
