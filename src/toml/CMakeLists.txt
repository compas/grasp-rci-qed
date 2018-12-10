cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)

add_library(toml STATIC
    toml.cc
    toml.F90
)
setup_fortran_modules(toml)
