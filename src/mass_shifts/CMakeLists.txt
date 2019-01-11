add_library(rciqed-mass_shifts STATIC
    mass_shifts.f90
    # Legacy code
    keint.f90 keint_I.f90
    vint.f90 vint_I.f90
    vinti.f90 vinti_I.f90
)
setup_fortran_modules(rciqed-mass_shifts)
target_link_libraries_Fortran(rciqed-mass_shifts PRIVATE
    mod 9290
    rciqed-grasp
)
