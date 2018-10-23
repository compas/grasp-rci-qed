add_library(rci_qed STATIC
    # Vacuum polarization related routines
    funk.f90 funk_I.f90
    funl.f90 funl_I.f90
    vac2.f90 vac2_I.f90
    vac4.f90 vac4_I.f90
    vacpol.f90 vacpol_I.f90
    # Self-energy related routines
    fzalf.f90 fzalf_I.f90
    klamaq.f90 klamaq_I.f90
    mohr.f90 mohr_I.f90
    qed.f90 qed_I.f90
    qed_slfen.f90 qed_slfen_I.f90
    ratden_I.f90
)
setup_fortran_modules(rci_qed)
target_link_libraries_Fortran(rci_qed PRIVATE mod 9290 mpi90)