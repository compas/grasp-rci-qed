# May need to set:
#
#    GFORTRAN_CONVERT_UNIT="big_endian"
#
using Test

length(ARGS) == 1 || error("Must pass compiled libgrasp-rci.so as first argument")
const libgrasp = abspath(first(ARGS))
isfile(libgrasp) || error("$(libgrasp) missing")

function grasp_initialize_constants()
    @debug "Starting ccall: grasp_initialize_constants()"
    status = ccall( (:grasp_initialize_constants, libgrasp), Cint, ())
    status == 0 || error("grasp_initialize_constants exited with an error")
    return
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

function grasp_orbital_grid(i::Integer)
    @debug "Starting ccall: grasp_ci_diracnuclear"
    r = ccall( (:grasp_orbital_grid, libgrasp), Cdouble, (Cint,), i)
    return Float64(r)
end

dc(i, j) = grasp_ci_diracnuclear(i,j) + grasp_ci_coulomb(i,j)

@testset "libgrasp-rci" begin
    @test_throws ErrorException grasp_load_isodata("isodata.Z60")
    @test grasp_initialize_constants() === nothing
    @test grasp_load_isodata("test/isodata") === nothing
    @test grasp_load_csl("test/test.c") === nothing
    # Note: requires isodata and CSL
    @test grasp_load_orbitals("test/test.w") === nothing
    @test grasp_init_rkco_gg() === nothing
    @test grasp_init_dcb() === nothing
    @test grasp_load_mixing("test/test.cm") === nothing
end

@show grasp_orbital_grid(1)
@show grasp_orbital_grid(2)
@show grasp_orbital_grid(3)
@show grasp_orbital_grid(4)

println("NCF = $(grasp_ncsfs())")

@show grasp_ci_coulomb(1, 1)
@show grasp_ci_diracnuclear(1, 1)

println("End of script.")
