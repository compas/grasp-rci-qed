add_library(rciqed-qed STATIC
    qed.f90
    # Legacy code
    ## Vacuum polarization related routines
    legacy/funk.f90 legacy/funk_I.f90
    legacy/funl.f90 legacy/funl_I.f90
    legacy/ncharg.f90 legacy/ncharg_I.f90
    legacy/vac2.f90 legacy/vac2_I.f90
    legacy/vac4.f90 legacy/vac4_I.f90
    legacy/vacpol.f90 legacy/vacpol_I.f90
    legacy/vpint.f90 legacy/vpint_I.f90
    legacy/vpintf.f90 legacy/vpintf_I.f90
    ## Self-energy related routines
    legacy/fzalf.f90 legacy/fzalf_I.f90
    legacy/klamaq.f90 legacy/klamaq_I.f90
    legacy/mohr.f90 legacy/mohr_I.f90
    legacy/qed_slfen.f90 legacy/qed_slfen_I.f90
    legacy/ratden.f90 legacy/ratden_I.f90
)
setup_fortran_modules(rciqed-qed)
target_link_libraries_Fortran(rciqed-qed PRIVATE mod 9290 rciqed-grasp)
