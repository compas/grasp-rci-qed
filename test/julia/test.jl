# This can be set with
#
#    process.env["GFORTRAN_CONVERT_UNIT"]="big_endian"
#
# in the developer tools of Atom, before starting the Julia process.
@assert ENV["GFORTRAN_CONVERT_UNIT"] === "big_endian"
using RCIWrapper
RCIWrapper.initialize!("test/isodata", "test/test.c", "test/test.w", "test/test.cm")

RCIWrapper.globals_orbitals()

op = RCIWrapper.qedse(2)
size(op)
RCIWrapper.materialize(op)

RCIWrapper.onescalar(1,1)

RCIWrapper.matrixelement(op, 1, 3)

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
