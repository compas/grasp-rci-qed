# Integration tests that run rci-qed in various configurations. The output .csum
# is checked byte-by-byte, against a
#
# TODO:
#
#   - A script to automatically update the reference output files.

# T1: Oxygen-like system, with 1s in core, 2s/2p-/2p in the valence.
add_test(NAME integration-T1 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/integration.sh "T1")
add_test(NAME integration-T2 COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/integration.sh "T2")

# Set up environment variables that the integration.sh script needs that are
# common to all the tests.
set_property(TEST integration-T1 integration-T2
    PROPERTY ENVIRONMENT
    "GRASP=${GRASP}"
    "GRASP_BUILD_BINDIR=${BUILD_BINDIR}"
    "GRASP_EXPORTHYDROGENIC=$<TARGET_FILE:rci-qed.exporthydrogenic>"
)
