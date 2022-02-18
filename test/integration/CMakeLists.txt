# Integration tests that run rci-qed in various configurations. The output .csum
# is checked byte-by-byte, against reference output.

set(INTEGRATIONSH "${CMAKE_CURRENT_SOURCE_DIR}/integration.sh")
set(INTEGRATION_ARGS "")
# In order to update the reference files, in case there are legitimate changes to
# to the output, the following line should be temporarily commented out and then
# one should run all the integration tests.
#set(INTEGRATION_ARGS "--update")
# To run just integration tests, the easiest way is to invoke CTest with:
#
#   ctest -R 'integration-*'
#

# T1 and T2 are both P-like W systems with the 3s(2)3p(3) configurations.
# T1 has no core, but T2 has n=1,2 as core.
add_test(NAME integration-T1 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T1")
add_test(NAME integration-T2 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T2")
# T3 is the same system, but has excitations up to 4s from 3s/3p. T3 is the same
# as T4, but runs with mpirun -n 5.
# TODO: there are minor differences in numbers between the T3 & T4 .csum files.
add_test(NAME integration-T3 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T3")
add_test(NAME integration-T4 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T4")
# T5 is same as T3, but limits the number of eigenvalues calculated in each block.
add_test(NAME integration-T5 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T5")

# T6 is neutral tungsten. T6b sets maximum self-energy orbital to n=2.
add_test(NAME integration-T6 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T6")
add_test(NAME integration-T6b COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T6b")
# T6c/d/e test the alternative self-energy operators instead
add_test(NAME integration-T6c COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T6c")
add_test(NAME integration-T6d COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T6d")
add_test(NAME integration-T6e COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T6e")

# T7 is neutral tungsten but no QED. T8 is the same, but runs with `mpirun -n 5`
add_test(NAME integration-T7 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T7")
add_test(NAME integration-T8 COMMAND ${INTEGRATIONSH} ${INTEGRATION_ARGS} "T8")

# Set up environment variables that the integration.sh script needs that are
# common to all the tests.
set_property(TEST
    integration-T1 integration-T2 integration-T3 integration-T4 integration-T5
    integration-T6 integration-T6b integration-T6c integration-T6d integration-T6e
    integration-T7 integration-T8
    PROPERTY ENVIRONMENT
    "GRASP=${GRASP}"
    "GRASP_BUILD_BINDIR=${BUILD_BINDIR}"
    "GRASP_EXPORTHYDROGENIC=$<TARGET_FILE:rci-qed.exporthydrogenic>"
)
