# To build the CMake-based bits you first need to set up the build directory
# (out of tree builds are preferred). For that run:
#
#    mkdir build/ && cd build/ && cmake ..
#
# And then under the build/ directory simply call
#
#    make
#
# which will compile and install all the libraries to lib/
#

cmake_minimum_required (VERSION 2.6)
project (grasp-rci-qed)
enable_language(Fortran)
enable_testing()

# rci-qed links agains the libraries in ${GRASP}/lib, so the user needs to
# specify the GRASP enviroment variable that would point to the GRASP root
# directory when first setting up the build directory.
#
# The value gets cached though, so it does not have to be specified every time
# CMake reconfigures itself.
#
# If the cached value and the environment variable values do not agree, the
# cached value gets used. To update the path to GRASP, you should remove it from
# the cache. Or just re-create the build directory.
if(DEFINED GRASP)
    if(NOT "$ENV{GRASP}" STREQUAL "")
        if(NOT "$ENV{GRASP}" STREQUAL "${GRASP}")
            message(WARNING "\$GRASP environment differs from cached value.\nUsing the cached value: ${GRASP}")
        endif()
    endif()
else()
    if("$ENV{GRASP}" STREQUAL "")
        message(FATAL_ERROR "\$GRASP environment variable not set.")
    endif()
    set(GRASP $ENV{GRASP} CACHE STRING "Path to GRASP")
endif()

# "Release" will be the default build type, which gives us optimization flags etc.
# The other relevant option would be "Debug", which disables optimizations and
# enables debugging symbols. The debug build can be enabled when setting up the
# build directory with CMake:
#
#     cmake -DCMAKE_BUILD_TYPE=Debug ..
#
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Release Debug."
        FORCE
    )
endif(NOT CMAKE_BUILD_TYPE)

# Find the LAPACK, BLAS and MPI libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)

# Command: target_link_libraries_Fortran(target mode libraries...)
#
# Similar to target_link_libraries(), but will also set up paths so that the
# compiler could fine the the Fortran .mod files from of the libraries. Unlike
# for the standard command, mode ( = PUBLIC, PRIVATE) is mandatory.
#
# Modified version of: https://stackoverflow.com/a/43918277/1601695
#
# Example:
#
#     target_link_libraries_Fortran(rcsfsplit PRIVATE mod 9290)
#
function(target_link_libraries_Fortran target mode)
    target_link_libraries(${target} ${mode} ${ARGN})
    foreach(lib IN LISTS ARGN)
        target_include_directories(${target} ${mode} $<TARGET_PROPERTY:${lib},Fortran_MODULE_DIRECTORY>)
    endforeach()
endfunction()

# We put the compiled binaries into the bin/ subdirectory of the build directory
# and libraries into the lib/ subdirectory.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin/")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib/")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib/")
# To install the binaries into the standard ${GRASP}/bin/ directory, you need to
# call `make install`.
set(CMAKE_INSTALL_PREFIX ${GRASP})

# Additional Fortran compiler flags.
#
# -fno-automatic: this was set in the original make_environment_gfortran_UBC file.
#
# Note: optimization should be enabled on the Release target automatically.
#
# If need be, you can also set up linker flags. E.g.:
#
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libgfortran")
#
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic")

# Find and set up the GRASP libraries. They are assumed to live in ${GRASP}/lib.
#
# Note: it does not appear to be possible to use target_link_libraries on imported
# targets. Hence the need to set IMPORTED_LINK_INTERFACE_LIBRARIES and LINK_FLAGS
# manually via set_property commands.
add_library(mod STATIC IMPORTED)
set_property(TARGET mod PROPERTY IMPORTED_LOCATION "${GRASP}/lib/libmod.a")
set_property(TARGET mod PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/libmod")

add_library(9290 STATIC IMPORTED)
set_property(TARGET 9290 PROPERTY IMPORTED_LOCATION "${GRASP}/lib/lib9290.a")
set_property(TARGET 9290 PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/lib9290")
set_property(TARGET 9290 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES
    mod ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_property(TARGET 9290 PROPERTY LINK_FLAGS ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})

add_library(mcp90 STATIC IMPORTED)
set_property(TARGET mcp90 PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/libmcp90")
set_property(TARGET mcp90 PROPERTY IMPORTED_LOCATION "${GRASP}/lib/libmcp90.a")
set_property(TARGET mcp90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES mod 9290)

add_library(rang90 STATIC IMPORTED)
set_property(TARGET rang90 PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/librang90")
set_property(TARGET rang90 PROPERTY IMPORTED_LOCATION "${GRASP}/lib/librang90.a")
set_property(TARGET rang90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES mod 9290 mcp90)

add_library(mpi90 STATIC IMPORTED)
set_property(TARGET mpi90 PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/mpi90")
set_property(TARGET mpi90 PROPERTY IMPORTED_LOCATION "${GRASP}/lib/libmpiu90.a")
set_property(TARGET mpi90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES ${MPI_Fortran_LIBRARIES})
set_target_properties(mpi90 PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)

add_library(dvd90 STATIC IMPORTED)
set_property(TARGET dvd90 PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/libdvd90")
set_property(TARGET dvd90 PROPERTY IMPORTED_LOCATION "${GRASP}/lib/libdvd90.a")
set_property(TARGET dvd90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES  mod)


message("Compiler flags etc. for this rci-qed build:")
message("* GRASP:                          ${GRASP}")
message("* CMAKE_BUILD_TYPE:               ${CMAKE_BUILD_TYPE}")
message("* CMAKE_Fortran_COMPILER:         ${CMAKE_Fortran_COMPILER}")
message("* CMAKE_Fortran_COMPILER_VERSION: ${CMAKE_Fortran_COMPILER_VERSION}")
message("* CMAKE_Fortran_FLAGS:            ${CMAKE_Fortran_FLAGS}")
message("* CMAKE_Fortran_FLAGS_RELEASE:    ${CMAKE_Fortran_FLAGS_RELEASE}")
message("* CMAKE_Fortran_FLAGS_DEBUG:      ${CMAKE_Fortran_FLAGS_DEBUG}")
message("* CMAKE_EXE_LINKER_FLAGS:         ${CMAKE_EXE_LINKER_FLAGS}")
message("* CMAKE_STATIC_LINKER_FLAGS:      ${CMAKE_STATIC_LINKER_FLAGS}")
message("* CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message("* CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message("* CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
message("* BLAS_LIBRARIES:                 ${BLAS_LIBRARIES}")
message("* BLAS_LINKER_FLAGS:              ${BLAS_LINKER_FLAGS}")
message("* LAPACK_LIBRARIES:               ${LAPACK_LIBRARIES}")
message("* LAPACK_LINKER_FLAGS:            ${LAPACK_LINKER_FLAGS}")
message("* MPI_Fortran_LIBRARIES:          ${MPI_Fortran_LIBRARIES}")
message("* MPI_Fortran_INCLUDE_PATH:       ${MPI_Fortran_INCLUDE_PATH}")
message("* MPI_Fortran_COMPILE_FLAGS:      ${MPI_Fortran_COMPILE_FLAGS}")
message("* MPI_Fortran_LINK_FLAGS:         ${MPI_Fortran_LINK_FLAGS}")

# Include the rci-qed build target
add_subdirectory("src/")

# Include unit and integration tests
add_subdirectory("test/")
