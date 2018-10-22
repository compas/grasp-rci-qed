add_test(NAME test-integration-rci
    COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/integration-rci.sh
)
set_property(TEST test-integration-rci PROPERTY ENVIRONMENT
    "GRASP=${GRASP}"
    "RCIQED=$<TARGET_FILE:rci-qed>"
)
