add_executable(test.matrixelements_hydrogenic
    matrixelements_hydrogenic.f90
)
target_link_libraries_Fortran(test.matrixelements_hydrogenic PUBLIC
    mod 9290
    rciqed-grasp rciqed-matrixelements rciqed-rci rciqed-breit rciqed-qed rciqed-mass_shifts
    rciqed-tools rciqed-test
)
add_test(matrixelements_hydrogenic_generate "../../tools/rci-qed.exporthydrogenic")
add_test(matrixelements_hydrogenic test.matrixelements_hydrogenic)
set_property(TEST matrixelements_hydrogenic PROPERTY ENVIRONMENT
    "RCIQED_TESTDATA=${TEST_SOURCE_DIR}/data"
)
set_property(TEST matrixelements_hydrogenic PROPERTY DEPENDS
    matrixelements_hydrogenic_generate
)

# rcisettings
add_executable(test.matrixelements_rcisettings_expfloat
    matrixelements_rcisettings_expfloat.f90
)
target_link_libraries_Fortran(test.matrixelements_rcisettings_expfloat PUBLIC
    rciqed-grasp rciqed-matrixelements
)
add_test(matrixelements_rcisettings_expfloat test.matrixelements_rcisettings_expfloat)
add_test(NAME matrixelements_rcisettings_expfloat-check
    COMMAND
    ${TEST_SOURCE_DIR}/lib/filematch.sh
    "expfloat.toml"
    ${CMAKE_CURRENT_SOURCE_DIR}/expfloat.toml.reference
)
set_tests_properties(matrixelements_rcisettings_expfloat
    PROPERTIES FIXTURES_SETUP rcisettings-expfloattoml)
set_tests_properties(matrixelements_rcisettings_expfloat-check
    PROPERTIES FIXTURES_REQUIRED rcisettings-expfloattoml)

# matrixelements_rcisettings
add_executable(test.matrixelements_rcisettings
    matrixelements_rcisettings.f90
)
target_link_libraries_Fortran(test.matrixelements_rcisettings PUBLIC
    mod
    rciqed-grasp rciqed-matrixelements
    rciqed-test
)
add_test(matrixelements_rcisettings test.matrixelements_rcisettings)
add_test(NAME matrixelements_rcisettings-check
    COMMAND
    ${TEST_SOURCE_DIR}/lib/filematch.sh
    "testjob.settings.toml"
    ${CMAKE_CURRENT_SOURCE_DIR}/testjob.settings.toml.reference
)
set_tests_properties(matrixelements_rcisettings
    PROPERTIES FIXTURES_SETUP rcisettings)
set_tests_properties(matrixelements_rcisettings-check
    PROPERTIES FIXTURES_REQUIRED rcisettings)
