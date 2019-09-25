# May need to set:
#
#    GFORTRAN_CONVERT_UNIT="big_endian"
#
using AtomicLevels
using Test
using Libdl

length(ARGS) == 1 || error("Must pass compiled libgrasp-rci.so as first argument")
const libgrasp = abspath(first(ARGS))
isfile(libgrasp) || error("$(libgrasp_path) missing")
const libgrasp_lib = Libdl.dlopen(libgrasp)

function grasp_initialize_constants()
    @debug "Starting ccall: grasp_initialize_constants()"
    sym = Libdl.dlsym(libgrasp_lib, :grasp_initialize_constants)
    status = ccall(sym, Cint, ())
    status == 0 || error("grasp_initialize_constants exited with an error")
    return
end

function grasp_libgraspci_j2max()
    @debug "Starting ccall: grasp_libgraspci_j2max"
    r = ccall( (:grasp_libgraspci_j2max, libgrasp), Cint, ())
    Int(r)
end

function grasp_error_string()
    @debug "Starting ccall: grasp_error_string"
    s = ccall( (:grasp_error_string, libgrasp), Cstring, ())
    unsafe_string(s)
end

for routine in [:grasp_load_isodata, :grasp_load_orbitals, :grasp_load_mixing, :grasp_load_csl]
    routine_literal = QuoteNode(routine)
    debugstr = Expr(:string, "Starting ccall: $(routine)('", :filename, "')")
    @show routine, routine_literal, debugstr
    expr = quote
        function $(routine)(filename)
            @debug $(debugstr)
            status = ccall( ($(routine_literal), libgrasp),
                Cint, (Cstring,),
                filename
            )
            if status != 0
                error(grasp_error_string())
            end
            return
        end
    end
    eval(expr)
end

function grasp_init_rkco_gg()
    @debug "Starting ccall: grasp_init_dcb"
    status = ccall( (:grasp_init_rkco_gg, libgrasp),
        Cint, (),
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

function grasp_init_dcb()
    @debug "Starting ccall: grasp_init_dcb"
    status = ccall( (:grasp_init_dcb, libgrasp),
        Cint, (),
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

function grasp_ncsfs()
    @debug "Starting ccall: grasp_ncsfs"
    ncsfs = ccall( (:grasp_ncsfs, libgrasp),
        Cint, (),
    )
    return Int(ncsfs)
end

function grasp_ci_init_qedvp()
    @debug "Starting ccall: grasp_ci_init_qedvp"
    ccall( (:grasp_ci_init_qedvp, libgrasp),
        Cvoid, (),
    )
    return nothing
end


function grasp_ci_coulomb(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_coulomb"
    r = ccall( (:grasp_ci_coulomb, libgrasp),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_ci_diracnuclear(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_diracnuclear"
    r = ccall( (:grasp_ci_diracnuclear, libgrasp),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_ci_qedvp(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_qedvp"
    r = ccall( (:grasp_ci_qedvp, libgrasp),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_orbital_grid(i::Integer)
    @debug "Starting ccall: grasp_orbital_grid"
    r = ccall( (:grasp_orbital_grid, libgrasp), Cdouble, (Cint,), i)
    return Float64(r)
end

struct OneScalarCache
    p :: Ptr{Cvoid}
end

function grasp_ci_onescalar(ic::Integer, ir::Integer)
    @debug "Starting ccall: grasp_ci_onescalar"
    r = ccall( (:grasp_ci_onescalar, libgrasp), Ptr{Cvoid}, (Cint,Cint), ic, ir)
    return OneScalarCache(r)
end

function grasp_ci_onescalar_show(c::OneScalarCache)
    @debug "Starting ccall: grasp_ci_onescalar_show"
    ccall( (:grasp_ci_onescalar_show, libgrasp), Cvoid, (Ptr{Cvoid},), c.p)
end

function grasp_ci_init(isodata, csl, orbitals, mixing)
    @debug "Starting ccall: grasp_ci_init" isodata csl orbitals mixing
    status = ccall( (:grasp_ci_init, libgrasp),
        Cint, (Cstring, Cstring, Cstring, Cstring),
        isodata, csl, orbitals, mixing
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

struct Matrix1PCache
    p :: Ptr{Cvoid}
end

function grasp_ci_qedse_matrix(setype::Integer)
    @debug "Starting ccall: grasp_ci_qedse_matrix"
    p = Ref{Ptr{Cvoid}}()
    status = ccall( (:grasp_ci_qedse_matrix, libgrasp),
        Cint, (Cint, Ptr{Ptr{Cvoid}}),
        setype, p
    )
    if status != 0
        error(grasp_error_string())
    end
    return Matrix1PCache(p[])
end

function grasp_ci_matrix1p(matrix1p::Matrix1PCache, ic::Integer, ir::Integer)
    @debug "Starting ccall: grasp_ci_matrix1p"
    osc = grasp_ci_onescalar(ic, ir)
    r = ccall( (:grasp_ci_matrixelement_1p_cached, libgrasp),
        Cdouble, (Ptr{Cvoid}, Ptr{Cvoid}),
        osc.p, matrix1p.p
    )
    Float64(r)
end

function grasp_libgraspci_matrix1p(matrix1p::Matrix1PCache)
    @debug "Starting ccall: grasp_ci_matrix1p"
    nw = grasp_orbitals_nw()
    elem(i, j) = ccall( (:grasp_libgraspci_matrix1p, libgrasp),
        Cdouble, (Ptr{Cvoid}, Cint, Cint),
        matrix1p.p, i, j
    )
    [elem(i, j) for i=1:nw, j=1:nw]
end

function grasp_orbitals_nw()
    @debug "Starting ccall: grasp_orbitals_nw"
    r = ccall( (:grasp_orbitals_nw, libgrasp), Cint, ())
    Int(r)
end

function grasp_orbitals()
    @debug "Starting ccall: grasp_orbitals"
    nw = grasp_orbitals_nw()
    np, nak = Vector{Cint}(undef, nw), Vector{Cint}(undef, nw)
    nw_ref = Ref{Cint}()
    ccall( (:grasp_orbitals, libgrasp),
        Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        nw_ref, np, nak
    )
    return map(args -> RelativisticOrbital(args...), zip(np, nak))
end

dc(i, j) = grasp_ci_diracnuclear(i, j) + grasp_ci_coulomb(i, j)

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

@show grasp_libgraspci_j2max()
grasp_ci_init("test/isodata", "test/test.c", "test/test.w", "test/test.cm")
@show grasp_libgraspci_j2max()

@show grasp_orbital_grid(1)
@show grasp_orbital_grid(2)
@show grasp_orbital_grid(3)
@show grasp_orbital_grid(4)

println("NCF = $(grasp_ncsfs())")

@show grasp_ci_coulomb(1, 1)
@show grasp_ci_diracnuclear(1, 1)

println("End of script.")

grasp_orbitals()
