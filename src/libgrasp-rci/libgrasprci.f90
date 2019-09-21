!! Contains the global state of libgrasprci
module libgrasprci

    ! Error handling
    character(:), allocatable :: lasterror
    character(:), allocatable, target :: lasterror_cstr

    logical :: constants_initialized = .false.
    logical :: isodata_loaded = .false.
    logical :: csl_loaded = .false.
    logical :: orbitals_loaded = .false.
    logical :: mixing_loaded = .false.

    integer :: j2max

end module libgrasprci
