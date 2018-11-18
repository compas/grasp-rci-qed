add_library(rci_breit STATIC
    bessel.f90 bessel_I.f90
    breid.f90 breid_I.f90
    brint1.f90 brint1_I.f90
    brint2.f90 brint2_I.f90
    brint3.f90 brint3_I.f90
    brint4.f90 brint4_I.f90
    brint5.f90 brint5_I.f90
    brint6.f90 brint6_I.f90
    brintf.f90 brintf_I.f90
    brra.f90 brra_I.f90
    cxk.f90 cxk_I.f90
    genintbreit1.f90 genintbreit1_I.f90
    genintbreit1wrap.f90 genintbreit1wrap_I.f90
    genintbreit2.f90 genintbreit2_I.f90
    genintbreit2wrap.f90 genintbreit2wrap_I.f90
    rkint.f90 rkint_I.f90
    skint.f90 skint_I.f90
    talk.f90 talk_I.f90
    triangbreit1.f90 triangbreit1_I.f90
    triangbreit2.f90 triangbreit2_I.f90
    zkf.f90 zkf_I.f90
)
setup_fortran_modules(rci_breit)
target_link_libraries_Fortran(rci_breit PRIVATE mod 9290 rang90)
