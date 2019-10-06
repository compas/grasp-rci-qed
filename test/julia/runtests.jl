using Test
using AtomicLevels
using RCIQED

const DATADIR = joinpath(@__DIR__, "data")

@testset "Hydrogenic" begin
    RCIQED.__reload__()
    RCIQED.initialize!(
        joinpath(DATADIR, "hydrogen/isodata"),
        joinpath(DATADIR, "hydrogen/hydrogen.c"),
        joinpath(DATADIR, "hydrogen/hydrogen.w"),
        joinpath(DATADIR, "hydrogen/hydrogen.cm"),
    )

    os, es = RCIQED.globals_orbitals()
    @test os == ros"1[s] 2[s-p] 3[s-d] 4[s-f]"
    @test es[1]  ≈ 0.5000066565 atol=1e-10 # 1s
    @test es[2]  ≈ 0.1250020801 atol=1e-10 # 2s
    @test es[3]  ≈ 0.1250020801 atol=1e-10 # 2p-
    @test es[4]  ≈ 0.1250004160 atol=1e-10 # 2p
    # TODO: These energies do not match.. is it because of screening?
    @test_broken es[10] ≈ 0.0555558020 atol=1e-9 # 3p
    @test_broken es[14] ≈ 0.0312501300 atol=1e-9 # 4d-

    @test RCIQED.globals(:orb_C, :NCF) == 16 # Number of CSFs
    @test RCIQED.globals(:prnt_C, :NVEC) == 16 # Number of ASFs
end

@testset "Be-like" begin
    RCIQED.__reload__()
    RCIQED.initialize!(
        joinpath(DATADIR, "belike/isodata"),
        joinpath(DATADIR, "belike/belike.c"),
        joinpath(DATADIR, "belike/belike.w"),
        joinpath(DATADIR, "belike/belike.cm"),
    )
    RCIQED.initialize_breit!()
    RCIQED.initialize_qedvp!()
    RCIQED.initialize_mass_shifts!()

    os, es = RCIQED.globals_orbitals()
    @test os == ros"1[s] 2[s-p] 3[s-d]"

    @test RCIQED.globals(:orb_C, :NCF) == 211 + 436 + 534 # Number of CSFs
    @test RCIQED.globals(:prnt_C, :NVEC) == 10 + 9 + 11 # Number of ASFs
end
