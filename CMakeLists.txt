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

include("cmake/fortran_modules.cmake")

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

# If true, looks for the Fortran .mod files in ${GRASP}/lib/${library}, instead of the
# default source directories (i.e. \${GRASP}/src/lib/\${library}).
#
# The original GRASP 2018 publication leaves the .mod files next to the .f90 source files
# in the source directories of the libraries, so we support that file organization by
# default.
#
# However, the new CMake-based build on the GRASP master will install the .mod files next
# to the archive files under ${GRASP}/lib. So, the user can pass -DGRASP_INSTALLED_MODULES=true
# to cmake to use that instead (mandatory when linking against a CMake-based GRASP build).
set(GRASP_INSTALLED_MODULES false CACHE BOOL
    "If true, looks for the Fortran .mod files in \${GRASP}/lib/\${library}, instead of the
    default source directories (i.e. \${GRASP}/src/lib/\${library})"
)

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
function(add_grasp_static_library target directory)
    add_library(${target} STATIC IMPORTED)
    set_property(TARGET ${target} PROPERTY IMPORTED_LOCATION "${GRASP}/lib/lib${target}.a")
    if(GRASP_INSTALLED_MODULES)
        set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/lib/${target}")
    else(GRASP_INSTALLED_MODULES)
        set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${GRASP}/src/lib/${directory}")
    endif(GRASP_INSTALLED_MODULES)
endfunction()
# libmod.a
add_grasp_static_library(mod "libmod")
get_property(FOO TARGET mod PROPERTY Fortran_MODULE_DIRECTORY)
# lib9290.a
add_grasp_static_library(9290 "lib9290")
set_property(TARGET 9290 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES
    mod ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
set_property(TARGET 9290 PROPERTY LINK_FLAGS ${BLAS_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
# libmcp90.a
add_grasp_static_library(mcp90 "libmcp90")
set_property(TARGET mcp90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES mod 9290)
# librang90.a
add_grasp_static_library(rang90 "librang90")
set_property(TARGET rang90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES mod 9290 mcp90)
# libmpiu90.a
add_grasp_static_library(mpiu90 "mpi90")
set_property(TARGET mpiu90 PROPERTY IMPORTED_LINK_INTERFACE_LIBRARIES ${MPI_Fortran_LIBRARIES})
set_target_properties(mpiu90 PROPERTIES
  COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
  LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
)
# libdvd90.a
add_grasp_static_library(dvd90 "libdvd90")
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
add_subdirectory("tools/")
add_subdirectory("test/")
