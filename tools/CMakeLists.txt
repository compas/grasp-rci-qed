# The libraries and executables from tools/ should not end up in the standard
# build/lib/ and build/bin/ directories.
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Abort compilation if implicit procedures or interfaces are being used
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-procedure -Werror=implicit-interface")

add_library(rciqed-tools STATIC
    lib/lib9290_hydrogenic.f90
    lib/lib9290_setup.f90
    lib/orbout_I.f90 lib/orbout.f90 # from rmcdhf90
)
setup_fortran_modules(rciqed-tools)
target_link_libraries_Fortran(rciqed-tools PUBLIC
    mod 9290 mpiu90
    rciqed-grasp
)

# exporthydrogenic
add_executable(rci-qed.exporthydrogenic
    exporthydrogenic.f90
)
target_link_libraries_Fortran(rci-qed.exporthydrogenic PUBLIC
    mod 9290 mpiu90
    rciqed-tools
)

# exportvp
add_executable(rci-qed.exportvp
    exportvp.f90
)
target_link_libraries_Fortran(rci-qed.exportvp PUBLIC
    mod 9290 mpiu90
    rciqed-rci rciqed-qed
)
