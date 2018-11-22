# Warn if implicit procedures or interfaces are being used
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wimplicit-interface -Wimplicit-procedure")

add_library(grasp STATIC
    kinds.f90
    lib9290_init.f90
    lib9290_csls.F90
    lib9290_files.f90
    lib9290.f90
    # Legacy routines. Should be only used internally in this library.
    legacy/getmixblock.f90 legacy/getmixblock_I.f90
)
setup_fortran_modules(grasp)
target_link_libraries_Fortran(grasp PRIVATE
    mod 9290
)
