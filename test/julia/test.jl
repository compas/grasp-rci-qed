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
using RCIWrapper
RCIWrapper.__reload__()
#RCIWrapper.initialize!("test/isodata", "test/test.c", "test/test.w", "test/test.cm")
RCIWrapper.initialize!("Z30.cc3-n3/isodata.Z30", "Z30.cc3-n3/F_like-cc-n3.c", "Z30.cc3-n3/Z30_F_like-cc-n3.w", "Z30.cc3-n3/Z30.cc-n3.dcb.cm")
RCIWrapper.initialize_breit!()
RCIWrapper.initialize_qedvp!()
RCIWrapper.initialize_mass_shifts!()

RCIWrapper.globals_orbitals()

RCIWrapper.globals(:orb_C, :NCF) # Number of CSFs
RCIWrapper.globals(:prnt_C, :NVEC) # Number of ASFs

# Compute the single-particle matrix elements of the QED SE operators
ops = [RCIWrapper.qedse(i) for i=1:4]


dpop = RCIWrapper.diracpot()
RCIWrapper.asfvalue(dpop, 1)
RCIWrapper.asfvalue(RCIWrapper.diracpot, 1)
RCIWrapper.asfvalue(RCIWrapper.coulomb, 1)
RCIWrapper.asfvalue(RCIWrapper.breit, 1)
RCIWrapper.asfvalue(RCIWrapper.nms, 1)
RCIWrapper.asfvalue(RCIWrapper.sms, 1)
RCIWrapper.asfvalue(RCIWrapper.qedvp, 1)

op = RCIWrapper.qedse(2)
size(op)
RCIWrapper.materialize(op)

ops = [
    RCIWrapper.diracpot, RCIWrapper.coulomb,
    RCIWrapper.breit,
    RCIWrapper.nms, RCIWrapper.sms,
    RCIWrapper.qedvp,
    (RCIWrapper.qedse(i) for i=1:4)...
]
RCIWrapper.asfvalues(ops)

RCIWrapper.onescalar(1,1)

RCIWrapper.matrixelement(op, 1, 1)

RCIWrapper.asfcoefficients()

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

RCIWrapper._reopen_libgrasp()
