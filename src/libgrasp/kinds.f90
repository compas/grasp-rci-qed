!> Small shorthand module providing the commonly used kinds.
!!
!! It contains the `real64` and `dp` parameters, which are kind parameters for
!! 64-bit floats. Variable declarations should therefore use `real(real64)` and
!! floating point literals should always have `_dp` appended (e.g. `0.23_dp`).
!!
module grasp_rciqed_kinds
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    integer, parameter :: dp = real64

end module grasp_rciqed_kinds
