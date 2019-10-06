module RCIQED
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
const libgrasp = Ref{String}("")
"Contains the `dlopen`ed pointer to the `libgrasp-rci` library."
const libgrasp_lib = Ref{Ptr{Nothing}}(C_NULL)

"""
    __reload__()

Convenience function to quickly reload the `libgrasp-rci` shared library. Note
that this will wipe out any global state you have.
"""
function __reload__()
    global libgrasp, libgrasp_lib
    # Check the LIBGRASPRCI environment variable
    libgrasp[] = if haskey(ENV, "LIBGRASPRCI")
        @info "Overriding library location with \$LIBGRASPRCI" ENV["LIBGRASPRCI"]
        normpath(abspath(expanduser(ENV["LIBGRASPRCI"])))
    else
        normpath(abspath(joinpath(@__DIR__, "..", "lib", "libgrasp-rci.so")))
    end
    isfile(libgrasp[]) || error("libgrasp-rci.so not found at $(libgrasp[])")
    # If libgrasp is already loaded, unload it
    if libgrasp_lib[] != C_NULL
        Libdl.dlclose(libgrasp_lib[])
        libgrasp_lib[] = C_NULL
    end
    # Open the shared library (again)
    libgrasp_lib[] = Libdl.dlopen(libgrasp[])
    return
end

__init__() = __reload__()


"""
    grasperror() -> String

Returns the error string corresponding to the last error condition in
`libgrasp-rci`.
"""
function grasperror()
    s = ccall(Libdl.dlsym(libgrasp_lib[], :libgrasprci_error_string), Cstring, ())
    unsafe_string(s)
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
        error(grasperror())
    end
    return
end

function initialize_breit!()
    status = ccall(Libdl.dlsym(libgrasp_lib[], :libgrasprci_initialize_breit),
        Cint, ()
    )
    if status != 0
        error(grasperror())
    end
end

function initialize_qedvp!()
    status = ccall(Libdl.dlsym(libgrasp_lib[], :libgrasprci_initialize_qedvp),
        Cint, ()
    )
    if status != 0
        error(grasperror())
    end
end

function initialize_mass_shifts!()
    status = ccall(Libdl.dlsym(libgrasp_lib[], :libgrasprci_initialize_ms),
        Cint, ()
    )
    if status != 0
        error(grasperror())
    end
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
    nw = globals(:orb_C, :NW)
    np, nak = Vector{Cint}(undef, nw), Vector{Cint}(undef, nw)
    e = Vector{Cdouble}(undef, nw)
    sym = Libdl.dlsym(libgrasp_lib[], :libgraspci_global_orbitals)
    ccall(sym, Cvoid, (Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), nw, np, nak, e)
    return map(args -> RelativisticOrbital(args...), zip(np, nak)), e
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
struct Coulomb <: Operator end
const coulomb = Coulomb()

"""
    struct Breit <: Operator

A singleton type representing the Breit operator.
"""
struct Breit <: Operator end
const breit = Breit()

"""
    struct NormalMassShift <: Operator

A singleton type representing the normal mass shift operator.
"""
struct NormalMassShift <: OneScalarOperator end
const nms = NormalMassShift()

"""
    struct SpecialMassShift <: Operator

A singleton type representing the special mass shift operator.
"""
struct SpecialMassShift <: Operator end
const sms = SpecialMassShift()

"""
    struct QEDVP <: Operator

A singleton type representing the QED vacuum polarization operator.
"""
struct QEDVP <: OneScalarOperator end
const qedvp = QEDVP()

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

function matrixelement(::NormalMassShift, osc::OneScalarCache)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_nms)
    r = ccall(sym, Cdouble, (Ptr{Cvoid},), osc.p)
    return Float64(r)
end

function matrixelement(::QEDVP, osc::OneScalarCache)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_qedvp)
    r = ccall(sym, Cdouble, (Ptr{Cvoid},), osc.p)
    return Float64(r)
end

function matrixelement(::Coulomb, ir::Integer, ic::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_coulomb)
    r = ccall(sym, Cdouble, (Cint, Cint), ir, ic)
    return Float64(r)
end

function matrixelement(::Breit, ir::Integer, ic::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_breit)
    r = ccall(sym, Cdouble, (Cint, Cint), ir, ic)
    return Float64(r)
end

function matrixelement(::SpecialMassShift, ir::Integer, ic::Integer)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_matrixelement_sms)
    r = ccall(sym, Cdouble, (Cint, Cint), ir, ic)
    return Float64(r)
end

#-------------------------------------------------------------------------------------------
# Calculations with Atomic State Functions (ASFs)
# -----------------------------------------------

function asfcoefficients()
    ncf, nvec = globals(:orb_C, :NCF), globals(:prnt_C, :NVEC)
    sym = Libdl.dlsym(libgrasp_lib[], :libgrasprci_asfcoefficient)
    asf(k, i)::Float64 = ccall(sym, Cdouble, (Cint, Cint), k, i)
    Float64[asf(k, i) for i = 1:ncf, k = 1:nvec]
end

function asfvalues(ops::Vector{T}) where {T <: Operator}
    ncf = globals(:orb_C, :NCF)
    asfs = asfcoefficients()
    asfvalues = zeros(Float64, length(ops), size(asfs, 2))
    is1ps = map(op -> isa(op, Matrix1PScalar), ops)
    cij = Vector{Float64}(undef, size(asfs, 2))
    for i = 1:ncf, j = 1:ncf
        osc = any(is1ps) ? onescalar(i, j) : nothing
        #cij[:] .= asfs[i,:] .* asfs[j,:]
        for k = 1:size(asfs, 2)
            cij[k] = asfs[i, k] * asfs[j, k]
        end
        for (q, op) in enumerate(ops)
            me = is1ps[q] ? matrixelement(op, osc) : matrixelement(op, i, j)
            #asfvalues[q,:] .+= cij  .* me
            for k = 1:size(asfs, 2)
                asfvalues[q, k] += cij[k] * me
            end
        end
    end
    return asfvalues
end

function asfvalue(op::Operator, k)
    ncf, asfs = globals(:orb_C, :NCF), asfcoefficients()
    sum(matrixelement(op, i, j) * asfs[i, k] * asfs[j, k] for i = 1:ncf, j=1:ncf)
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
