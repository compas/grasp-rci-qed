module RCIWrapper
using AtomicLevels
using Libdl

#-------------------------------------------------------------------------------------------
# Managing the libgrasp-rci dynamic library
# -----------------------------------------
# We use dlopen to open library before so that we could also reload the library on demand,
# which is handy during development, as we need to recompile the library often.
#
# `libgrasp` contains the filesystem path to the library and `libgrasp_lib` is the pointer

"Contains the filesystem path to the `libgrasp-rci` library."
const libgrasp = Ref{String}()
"Contains the `dlopen`ed pointer to the `libgrasp-rci` library."
const libgrasp_lib = Ref{Ptr{Nothing}}(0)

function __init__()
    haskey(ENV, "LIBGRASPRCI") || error("\$LIBGRASPRCI unset")
    libgrasp_path = normpath(abspath(expanduser(ENV["LIBGRASPRCI"])))
    isfile(libgrasp_path) || error("$(libgrasp_path) missing")
    libgrasp[] = libgrasp_path
    _open_libgrasp()
end

function _open_libgrasp()
    libgrasp_lib[] == Ptr{Nothing}(0) || error("libgrasp is already open")
    libgrasp_lib[] = Libdl.dlopen(libgrasp[])

    # Initialize GRASP constants
    ccall(Libdl.dlsym(libgrasp_lib[], :libgrasprci_initialize_constants), Cvoid, ())
    return
end

"""
Convenience function to quickly reload the `libgrasp-rci` library.
"""
function _reopen_libgrasp()
    Libdl.dlclose(libgrasp_lib[])
    libgrasp_lib[] = 0
    _open_libgrasp()
end

#-------------------------------------------------------------------------------------------
# Managing the libgrasp-rci global state
# --------------------------------------

"""
Load all the relevant GRASP input files.
"""
function initialize!(isodata, csl, orbitals, mixing)
    @debug "Starting ccall: grasp_ci_init" isodata csl orbitals mixing
    isfile(isodata) || error("$(isodata) not a file")
    isfile(csl) || error("$(csl) not a file")
    isfile(orbitals) || error("$(orbitals) not a file")
    isfile(mixing) || error("$(mixing) not a file")
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_initalize)
    status = ccall(sym,
        Cint, (Cstring, Cstring, Cstring, Cstring),
        isodata, csl, orbitals, mixing
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

function globals(mod::Symbol, variable::Symbol)
    if mod === :libgrasprci
        if variable === :j2max
            sym = Libdl.dlsym(libgrasp_lib[], :libgraspci_global_j2max)
            Int(ccall(sym, Cint, ()))
        else
            error("Unsupported variable in orb_C: $variable")
        end
    elseif mod === :orb_C
        if variable === :NW
            sym = Libdl.dlsym(libgrasp_lib[], :libgraspci_global_orb_nw)
            Int(ccall(sym, Cint, ()))
        elseif variable === :NCF
            sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_global_orb_ncf)
            Int(ccall(sym, Cint, ()))
        else
            error("Unsupported variable in orb_C: $variable")
        end
    elseif mod === :prnt_C
        if variable === :NVEC
            sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_global_prnt_nvec)
            Int(ccall(sym, Cint, ()))
        else
            error("Unsupported variable in prnt_C: $variable")
        end
    else
        error("Unsupported module $mod")
    end
end

"""
    globals_orbitals()

Returns a vector of `RelativisticOrbital` of the orbitals stored in the GRASP global state
(i.e. in the `orb_C` module).
"""
function globals_orbitals()
    nw = RCIWrapper.globals(:orb_C, :NW)
    np, nak = Vector{Cint}(undef, nw), Vector{Cint}(undef, nw)
    sym = Libdl.dlsym(libgrasp_lib[], :libgraspci_global_orbitals)
    ccall(sym, Cvoid, (Cint, Ptr{Cint}, Ptr{Cint}), nw, np, nak)
    return map(args -> RelativisticOrbital(args...), zip(np, nak))
end

function asfcoefficients()
    ncf, nvec = globals(:orb_C, :NCF), globals(:prnt_C, :NVEC)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_asfcoefficient)
    asf(k, i)::Float64 = ccall(sym, Cdouble, (Cint, Cint), k, i)
    Float64[asf(k, i) for i = 1:ncf, k = 1:nvec]
end


#-------------------------------------------------------------------------------------------
# 1 particle scalar operators (generic API)
# -----------------------------------------

abstract type Operator end
abstract type OneScalarOperator <: Operator end

"""
    mutable struct Matrix1PScalar

Stores a reference to an object of the Fortran `matrix1pscalar` type.
"""
mutable struct Matrix1PScalar <: OneScalarOperator
    p :: Ptr{Cvoid}

    function Matrix1PScalar(p::Ptr)
        x = new(p)
        finalizer(Matrix1PScalar_finalizer, x)
        return x
    end
end

function Matrix1PScalar_finalizer(m1ps::Matrix1PScalar)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar_delete)
    ccall(sym, Cvoid, (Ptr{Cvoid},), m1ps.p)
    return
end

function Base.size(m1ps::Matrix1PScalar)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar_n)
    n = ccall(sym, Cint, (Ptr{Cvoid},), m1ps.p)
    (Int(n), Int(n))
end

function Base.size(m1ps::Matrix1PScalar, k::Int)
    k <= 0 && error("dimension $(k) out of range")
    k > 2 && return 1
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar_n)
    n = ccall(sym, Cint, (Ptr{Cvoid},), m1ps.p)
    Int(n)
end

function Base.copy(m1ps::Matrix1PScalar)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar_copy)
    p = ccall(sym, Ptr{Cvoid}, (Ptr{Cvoid},), m1ps.p)
    Matrix1PScalar(p)
end

function materialize(m1ps::Matrix1PScalar)
    nw = size(m1ps, 1)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar)
    [
        ccall(sym, Cdouble, (Ptr{Cvoid}, Cint, Cint), m1ps.p, i, j)
        for i=1:nw, j=1:nw
    ]
end

function disable_orbital!(m1ps::Matrix1PScalar, i::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrix1pscalar_disable)
    ccall(sym, Cvoid, (Ptr{Cvoid}, Cint), m1ps.p, i)
    return
end

"""
    struct DiracPotential <: Operator

A singleton type representing the Dirac + central potential operator.
"""
struct DiracPotential <: OneScalarOperator end
const diracpot = DiracPotential()

"""
    struct Coulomb <: Operator

A singleton type representing the Coulomb operator.
"""
struct Coulomb <: OneScalarOperator end
const coulomb = Coulomb()

#-------------------------------------------------------------------------------------------
# libgrasp-rci QED operators
# --------------------------

function diracpot()
    p = Ref{Ptr{Cvoid}}()
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_diracpot_matrix1pscalar)
    status = ccall(sym, Cint, (Ptr{Ptr{Cvoid}},), p)
    if status != 0
        error(grasp_error_string())
    end
    return Matrix1PScalar(p[])
end

function qedse(setype::Integer)
    @debug "Starting ccall: grasp_ci_qedse_matrix"
    p = Ref{Ptr{Cvoid}}()
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_qedse_matrix1pscalar)
    status = ccall(sym, Cint, (Cint, Ptr{Ptr{Cvoid}}), setype, p)
    if status != 0
        error(grasp_error_string())
    end
    return Matrix1PScalar(p[])
end

#-------------------------------------------------------------------------------------------
# Matrix element calculation
# --------------------------

"""
    mutable struct OneScalarCache

Stores a reference to an object of the Fortran `onescalar_cache` type.
"""
mutable struct OneScalarCache
    p :: Ptr{Cvoid}

    function OneScalarCache(p::Ptr)
        x = new(p)
        finalizer(OneScalarCache_finalizer, x)
        return x
    end
end

function OneScalarCache_finalizer(osc::OneScalarCache)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_onescalar_delete)
    ccall(sym, Cvoid, (Ptr{Cvoid},), osc.p)
    return
end

function onescalar(ir::Integer, ic::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_onescalar)
    r = ccall(sym, Ptr{Cvoid}, (Cint, Cint), ir, ic)
    return OneScalarCache(r)
end

matrixelement(operator::OneScalarOperator, ir::Integer, ic::Integer) = matrixelement(operator, onescalar(ic, ir))

function matrixelement(operator::Matrix1PScalar, osc::OneScalarCache)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_1p)
    r = ccall(sym, Cdouble, (Ptr{Cvoid}, Ptr{Cvoid}), osc.p, operator.p)
    return Float64(r)
end

function matrixelement(::DiracPotential, osc::OneScalarCache)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_diracpot)
    r = ccall(sym, Cdouble, (Ptr{Cvoid},), osc.p)
    return Float64(r)
end

function matrixelement(::Coulomb, ir::Integer, ic::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_coulomb)
    r = ccall(sym, Cdouble, (Cint, Cint), ir, ic)
    return Float64(r)
end

#
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------

function grasp_ci_onescalar_show(c)
    @debug "Starting ccall: grasp_ci_onescalar_show"
    ccall( (:grasp_ci_onescalar_show, libgrasp[]), Cvoid, (Ptr{Cvoid},), c.p)
end

function grasp_error_string()
    @debug "Starting ccall: grasp_error_string"
    s = ccall( (:grasp_error_string, libgrasp[]), Cstring, ())
    unsafe_string(s)
end

for routine in [:grasp_load_isodata, :grasp_load_orbitals, :grasp_load_mixing, :grasp_load_csl]
    routine_literal = QuoteNode(routine)
    debugstr = Expr(:string, "Starting ccall: $(routine)('", :filename, "')")
    #@show routine, routine_literal, debugstr
    expr = quote
        function $(routine)(filename)
            @debug $(debugstr)
            status = ccall( ($(routine_literal), libgrasp[]),
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
    status = ccall( (:grasp_init_rkco_gg, libgrasp[]),
        Cint, (),
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

function grasp_init_dcb()
    @debug "Starting ccall: grasp_init_dcb"
    status = ccall( (:grasp_init_dcb, libgrasp[]),
        Cint, (),
    )
    if status != 0
        error(grasp_error_string())
    end
    return
end

function grasp_ci_init_qedvp()
    @debug "Starting ccall: grasp_ci_init_qedvp"
    ccall( (:grasp_ci_init_qedvp, libgrasp[]),
        Cvoid, (),
    )
    return nothing
end


function grasp_ci_coulomb(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_coulomb"
    r = ccall( (:grasp_ci_coulomb, libgrasp[]),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_ci_diracnuclear(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_diracnuclear"
    r = ccall( (:grasp_ci_diracnuclear, libgrasp[]),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_ci_qedvp(i::Integer, j::Integer)
    @debug "Starting ccall: grasp_ci_qedvp"
    r = ccall( (:grasp_ci_qedvp, libgrasp[]),
        Cdouble, (Cint, Cint),
        i, j
    )
    return Float64(r)
end

function grasp_orbital_grid(i::Integer)
    @debug "Starting ccall: grasp_orbital_grid"
    r = ccall( (:grasp_orbital_grid, libgrasp[]), Cdouble, (Cint,), i)
    return Float64(r)
end

end
