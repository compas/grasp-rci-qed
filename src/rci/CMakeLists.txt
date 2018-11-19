add_library(rci_rci
    auxblk.f90 auxblk_I.f90
    dnicmv.f90 dnicmv_I.f90
    dspevx_I.f90 # this interface appears to be for a LAPACK function
    engout.f90 engout_I.f90
    #gdsummpi_I.f90
    genintrk.f90 genintrk_I.f90
    genintrkwrap.f90 genintrkwrap_I.f90
    genmat.f90 genmat_I.f90
    genmat2.f90 genmat2_I.f90
    getcid.f90 getcid_I.f90
    hmout.f90 hmout_I.f90
    hovlap.f90 hovlap_I.f90
    iabint.f90 iabint_I.f90
    indtpi.f90 indtpi_I.f90
    iniestdm.f90 iniestdm_I.f90
    iniestsd.f90 iniestsd_I.f90
    lodmixmpi.f90 lodmixmpi_I.f90
    lodres.f90 lodres_I.f90
    maneig.f90 maneig_I.f90
    matrix.f90 matrix_I.f90
    rkintc.f90 rkintc_I.f90
    setdbg.f90 setdbg_I.f90
    setham_gg.f90 setham_gg_I.f90
    setmixmpi.f90 setmixmpi_I.f90
    setres.f90 setres_I.f90
    setsum.f90 setsum_I.f90
    spodmv.f90 spodmv_I.f90
    strsum.f90 strsum_I.f90
    triangrk.f90 triangrk_I.f90
    wghtd5.f90 wghtd5_I.f90
    # Added routines
    rci_utils.f90
)
setup_fortran_modules(rci_rci)
target_link_libraries_Fortran(rci_rci PRIVATE
    mod 9290 rang90 mpi90 dvd90
    rci_breit rci_qed rci_mass_shifts
)
