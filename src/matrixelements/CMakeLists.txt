# Warn if implicit interfaces are boung used.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_library(rciqed-matrixelements STATIC
    cimatrixelements.f90
    rcisettings.f90
    # Legacy code
    iabint.f90 iabint_I.f90
    rkintc.f90 rkintc_I.f90
)
setup_fortran_modules(rciqed-matrixelements)
target_link_libraries_Fortran(rciqed-matrixelements PRIVATE
    mod 9290 rang90
    rciqed-grasp rciqed-breit rciqed-qed rciqed-mass_shifts rciqed-toml
)
