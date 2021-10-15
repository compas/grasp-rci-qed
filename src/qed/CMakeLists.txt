add_library(rciqed-qed STATIC
    qed.f90
    pyykkoe.f90
    flambaum.f90
    qedmod.f90
    vp.f90

    # Legacy code
    ## Vacuum polarization related routines
    legacy/ncharg.f90 legacy/ncharg_I.f90
    legacy/vac2.f90 legacy/vac2_I.f90
    legacy/vac4.f90 legacy/vac4_I.f90
    legacy/vacpol.f90 legacy/vacpol_I.f90
    ## Self-energy related routines
    legacy/fzalf.f90 legacy/fzalf_I.f90
    legacy/klamaq.f90 legacy/klamaq_I.f90
    legacy/mohr.f90 legacy/mohr_I.f90
    legacy/qed_slfen.f90 legacy/qed_slfen_I.f90
    legacy/ratden.f90 legacy/ratden_I.f90

    # QEDMOD routines
    qedmod/dirac.f qedmod/exit.f qedmod/fse_dat.f
    qedmod/grid.f90 qedmod/grid_args.f qedmod/init_se.f qedmod/integr.f
    qedmod/nucl.f qedmod/populate_hydrogenics.f90 qedmod/pot_wav.f
    qedmod/read_func.f qedmod/se_pot.f qedmod/uehling.f
    qedmod/wk.f qedmod/write_func.f qedmod/wfunc_interpolate.f
)
setup_fortran_modules(rciqed-qed)
target_link_libraries_Fortran(rciqed-qed PRIVATE mod 9290 rciqed-grasp)
target_include_directories(rciqed-qed PRIVATE "qedmod/")
