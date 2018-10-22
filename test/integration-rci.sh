#!/bin/bash
export SCRIP_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -z ${GRASP+x} ]; then >&2 echo "ERROR: \$GRASP variable is unset."; exit 1; fi
if [ -z ${RCIQED+x} ]; then >&2 echo "ERROR: \$RCIQED variable is unset."; exit 1; fi

# We'll need to call the rwnfestimate binary from GRASP
RWFNESTIMATE=${GRASP}/bin/rwfnestimate
if ! [ -f "${RWFNESTIMATE}" ]; then
	>&2 echo "ERROR: rwnfestimate binary missing."
	>&2 echo "  searched at ${RWFNESTIMATE}"
	exit 1
fi

if ! [ -f "${RCIQED}" ]; then
	>&2 echo "ERROR: rci-qed binary missing."
	>&2 echo "  searched at ${RCIQED}"
	exit 1
fi

# The test will run in the integration-rci directory, under build/test.
if [ -e "integration-rci" ]; then
	if ! [ -d "integration-rci" ]; then
		>&2 echo "integration-rci exists at ${PWD}, but not directory. Bailing."
		exit 2
	fi

	>&2 echo "WARNING: removing integration-rci at ${PWD}"
	rm -Rv integration-rci
fi
mkdir "integration-rci" && cd "integration-rci" || {
	>&2 echo "ERROR: failed to create integration-rci/ at ${PWD}"
	exit 2
}
# The rest of this script has the working directory set to integration-rci/
>&2 echo "INFO: Working directory $PWD"
>&2 echo "INFO: Test script in $SCRIP_DIR"

# Set up the input files for rwnfestimate
cp ${SCRIP_DIR}/nitrogen.isodata isodata || exit 1
cp ${SCRIP_DIR}/nitrogen.c rcsf.inp || exit 1
sed 's/[ \t]*#.*//' <<-EOF > rwfnestimate.input
	y                  # Default settings?
	2                  # Wavefunction type? 2 = Thomas-Fermi
	*                  # List of subshells? All.
EOF
${RWFNESTIMATE} < rwfnestimate.input || {
	>&2 echo "FATAL ERROR: rwfnestimate failed"
	exit 1
}


# Set up input files for RCI
statename="nitrogen"
cp rcsf.inp ${statename}.c
cp rwfn.inp ${statename}.w
sed 's/[ \t]*#.*//' <<-EOF > rci.input
	y                  # Default settings?
	${statename}       # Name of state
	y                  # Include Breit?
	y                  # Modify frequency?
	1e-6               # Frequency scaling factor
	y                  # Include vaccuum polarization?
	y                  # Include NMS?
	y                  # Include SMS?
	y                  # Estimate self-energy?
	8                  # Max n for self-energy?
	1                  # J blocks 1-3
	1-3
	1
EOF

# Using a temporary directory in order to not run into the bug where GRASP limits
# the length of the MPI_TMP variable.
export MPI_TMP=`mktemp -d`
echo "MPI_TMP=${MPI_TMP}"
mpirun -n 1 ${RCIQED} < rci.input || {
	>&2 echo "FATAL ERROR: rci-qed failed with $?"
	>&2 echo "INFO: Cleanup:" && rm -Rv "${MPI_TMP}"
	exit 1
}
>&2 echo "INFO: Cleanup:" && rm -Rv "${MPI_TMP}"

# If we get this far, we're happy.
>&2 echo "${RCIQED} ran successfully!"
exit 0
