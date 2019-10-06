# The test wavefunctions are stored in big-endian encoding. The shared library
# can be forced into big-endian mode by setting
#
#    export GFORTRAN_CONVERT_UNIT="big_endian"
#
# before starting Julia on the command line, or with
#
#    process.env["GFORTRAN_CONVERT_UNIT"]="big_endian"
#
# in the developer tools of Atom, and then restarting the Julia process.
#
using Test
@test ENV["GFORTRAN_CONVERT_UNIT"] === "big_endian"
using RCIQED
RCIQED.__reload__()

#RCIQED.initialize!("test/isodata", "test/test.c", "test/test.w", "test/test.cm")
#RCIQED.initialize!("Z30.cc3-n3/isodata.Z30", "Z30.cc3-n3/F_like-cc-n3.c", "Z30.cc3-n3/Z30_F_like-cc-n3.w", "Z30.cc3-n3/Z30.cc-n3.dcb.cm")
RCIQED.initialize!(
    "Z30.cc-n7.dcb/isodata.Z30", "Z30.cc-n7.dcb/F_like-cc-n7.c",
    "Z30.cc-n7.dcb/Z30_F_like-cc-n7.w", "Z30.cc-n7.dcb/Z30.cc-n7.dcb.cm"
)
RCIQED.initialize_breit!()
RCIQED.initialize_qedvp!()
RCIQED.initialize_mass_shifts!()

RCIQED.globals_orbitals()

RCIQED.globals(:orb_C, :NCF) # Number of CSFs
RCIQED.globals(:prnt_C, :NVEC) # Number of ASFs

# Compute the single-particle matrix elements of the QED SE operators
ops = [RCIQED.qedse(i) for i=1:4]


dpop = RCIQED.diracpot()
RCIQED.asfvalue(dpop, 1)
RCIQED.asfvalue(RCIQED.diracpot, 1)
RCIQED.asfvalue(RCIQED.coulomb, 1)
RCIQED.asfvalue(RCIQED.breit, 1)
RCIQED.asfvalue(RCIQED.nms, 1)
RCIQED.asfvalue(RCIQED.sms, 1)
RCIQED.asfvalue(RCIQED.qedvp, 1)

op = RCIQED.qedse(2)
size(op)
RCIQED.materialize(op)

ops = [
    RCIQED.diracpot, RCIQED.coulomb,
    RCIQED.breit,
    RCIQED.nms, RCIQED.sms,
    RCIQED.qedvp,
    (RCIQED.qedse(i) for i=1:4)...
]
RCIQED.asfvalues(ops)

RCIQED.onescalar(1,1)

RCIQED.matrixelement(ops[2], 1, 1)

RCIQED.asfcoefficients()

# @testset "libgrasp-rci" begin
#     @test_throws ErrorException grasp_load_isodata("isodata.Z60")
#     @test grasp_initialize_constants() === nothing
#     @test grasp_load_isodata("test/isodata") === nothing
#     @test grasp_load_csl("test/test.c") === nothing
#     # Note: requires isodata and CSL
#     @test grasp_load_orbitals("test/test.w") === nothing
#     @test grasp_init_rkco_gg() === nothing
#     @test grasp_init_dcb() === nothing
#     @test grasp_load_mixing("test/test.cm") === nothing
# end

#@show grasp_libgraspci_j2max()
#@show grasp_libgraspci_j2max()

# @show grasp_orbital_grid(1)
# @show grasp_orbital_grid(2)
# @show grasp_orbital_grid(3)
# @show grasp_orbital_grid(4)
#
# println("NCF = $(grasp_ncsfs())")
#
# @show grasp_ci_coulomb(1, 1)
# @show grasp_ci_diracnuclear(1, 1)
#
# println("End of script.")
#
# grasp_orbitals()



# @show grasp_orbitals()
# m = grasp_ci_qedse_matrix(2)

RCIQED._reopen_libgrasp()
