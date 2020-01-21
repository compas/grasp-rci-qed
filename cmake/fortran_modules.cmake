# We need special functions to handle linking Fortran modules between libraries
# etc. The Fortran_MODULE_DIRECTORY_root variable is the directory where all the
# .mod files get written to. It is set be the modules/ subdirectory of the build
# directory.
#
# For every library, the modules get stored in
# ${Fortran_MODULE_DIRECTORY_root}/<library_name>/ so the modules from different
# libraries are separated from each other.
set(Fortran_MODULE_DIRECTORY_root ${CMAKE_CURRENT_BINARY_DIR}/modules)

# Command: setup_fortran_modules(target)
#
# Needs to be called on all libraries that provide modules. It set the
# Fortran_MODULE_DIRECTORY variable for the target, which is then used by
# target_link_libraries_Fortran to set up the appropriate include directories.
#
# Example:
#
#     setup_fortran_modules(9290)
#
function(setup_fortran_modules target)
    set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
endfunction()

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
