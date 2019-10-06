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
@assert ENV["GFORTRAN_CONVERT_UNIT"] === "big_endian"
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

# Operators
ops = [
    RCIQED.diracpot, RCIQED.coulomb,
    RCIQED.breit,
    RCIQED.nms, RCIQED.sms,
    RCIQED.qedvp,
    (RCIQED.qedse(i) for i=1:4)...
]
RCIQED.asfvalues(ops)


# ==============================================================================

RCIQED.matrixelement(ops[3], 1, 1)

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

RCIQED.onescalar(1,1)

RCIQED.asfcoefficients()
