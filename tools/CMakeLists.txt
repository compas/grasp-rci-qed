# The libraries and executables from tools/ should not end up in the standard
# build/lib/ and build/bin/ directories.
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)
unset(CMAKE_LIBRARY_OUTPUT_DIRECTORY)
unset(CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

add_library(rci-qed.tools STATIC
    lib/lib9290_hydrogenic.f90
    lib/lib9290_setup.f90
    lib/orbout_I.f90 lib/orbout.f90
)
setup_fortran_modules(rci-qed.tools)
target_link_libraries_Fortran(rci-qed.tools PUBLIC
    mod 9290 mpi90
    grasp
)

# exporthydrogenic
add_executable(rci-qed.exporthydrogenic
    exporthydrogenic.f90
)
target_link_libraries_Fortran(rci-qed.exporthydrogenic PUBLIC
    mod 9290 mpi90
    rci-qed.tools
)
#install(TARGETS rci-qed.exporthydrogenic DESTINATION bin/)
