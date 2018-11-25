#!/bin/bash
export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ -z ${GRASP+x} ]; then >&2 echo "ERROR: \$GRASP variable is unset."; exit 1; fi
if [ -z ${GRASP_BUILD_BINDIR+x} ]; then >&2 echo "ERROR: \$GRASP_BUILD_BINDIR variable is unset."; exit 1; fi

function find-grasp-binary {
	if [ "$#" -ne 1 ]; then
		>&2 echo "ERROR[find-grasp-binary]: unable to determine \$exe, bad arguments ($# '$@')"
		exit 1
	fi
	exe=$1
	guesses=(
		"${GRASP_BUILD_BINDIR}/${exe}"
		"${GRASP}/bin/${exe}"
	)
	for path in ${guesses[*]}; do
		if [ -f "$path" ]; then
			echo "$path"
			return
		fi
	done

	>&2 echo "ERROR[find-grasp-binary]: unable to find ${exe}"
	>&2 echo "  searched at:"
	for path in ${guesses[*]}; do
		>&2 echo "    ${path}"
	done
	exit 1
}

# We'll need to call the rwnfestimate binary from GRASP
RWFNESTIMATE=`find-grasp-binary rwfnestimate` || exit
RCIQED=`find-grasp-binary rci-qed` || exit

echo "RWFNESTIMATE: ${RWFNESTIMATE}"
echo "RCIQED:       ${RCIQED}"

# Get the configuration -- the first argument, and check that input files exist
CONF=$1
if ! [ -f "${SCRIPT_DIR}/${CONF}.isodata" ]; then >&2 echo "${SCRIPT_DIR}/${CONF}.isodata missing"; exit 1; fi
if ! [ -f "${SCRIPT_DIR}/${CONF}.c" ]; then >&2 echo "${SCRIPT_DIR}/${CONF}.c missing"; exit 1; fi

# The test will run in the integration-rci directory, under build/test.
TEST_DIRECTORY="integration-rci-${CONF}"
if [ -e "${TEST_DIRECTORY}" ]; then
	if ! [ -d "${TEST_DIRECTORY}" ]; then
		>&2 echo "${TEST_DIRECTORY} exists at ${PWD}, but not directory. Bailing."
		exit 2
	fi

	>&2 echo "WARNING: removing ${TEST_DIRECTORY} at ${PWD}"
	rm -Rv "${TEST_DIRECTORY}"
fi
mkdir "${TEST_DIRECTORY}" && cd "${TEST_DIRECTORY}" || {
	>&2 echo "ERROR: failed to create ${TEST_DIRECTORY}/ at ${PWD}"
	exit 2
}
# The rest of this script has the working directory set to integration-rci/
>&2 echo "INFO: Working directory $PWD"
>&2 echo "INFO: Test script in $SCRIPT_DIR"

# Set up the input files for rwnfestimate
cp ${SCRIPT_DIR}/${CONF}.isodata isodata || exit 1
cp ${SCRIPT_DIR}/${CONF}.c rcsf.inp || exit 1
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
cp rcsf.inp ${CONF}.c
cp rwfn.inp ${CONF}.w
sed 's/[ \t]*#.*//' <<-EOF > rci.input
	y                  # Default settings?
	${CONF}            # Name of state
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
