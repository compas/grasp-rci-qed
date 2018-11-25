#!/bin/bash
export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/../lib/testing.sh" || { echo "ERROR: Unable to load lib/testing.sh"; exit 1; }

if [ -z ${GRASP+x} ]; then >&2 echo "ERROR: \$GRASP variable is unset."; exit 1; fi
if [ -z ${GRASP_BUILD_BINDIR+x} ]; then >&2 echo "ERROR: \$GRASP_BUILD_BINDIR variable is unset."; exit 1; fi
if [ -z ${GRASP_EXPORTHYDROGENIC+x} ]; then >&2 echo "ERROR: \$GRASP_EXPORTHYDROGENIC variable is unset."; exit 1; fi
if ! [ -f ${GRASP_EXPORTHYDROGENIC} ]; then >&2 echo "ERROR: \$GRASP_EXPORTHYDROGENIC not pointing to a file."; exit 1; fi

# We'll need to locate the rci-qed binary
RCIQED=`find-grasp-binary rci-qed` || exit
echo "INFO: RCIQED=${RCIQED}"

# Get the configuration -- the first argument, and check that input files exist
TEST_CONF=$1
TEST_DIR="${SCRIPT_DIR}/$1"
check_directories "${TEST_DIR}"
check_files "${TEST_DIR}/isodata" "${TEST_DIR}/rci.input" "${TEST_DIR}/${TEST_CONF}.c"

# Test-specific configuration
if [ -f "${TEST_DIR}/settings.sh" ]; then
	source "${TEST_DIR}/settings.sh"
fi
GRASP_MPI_NPROCS=${GRASP_MPI_NPROCS:-1}
echo "GRASP_MPI_NPROCS=${GRASP_MPI_NPROCS}"

# The test will run in the integration-rci directory, under build/test.
TEST_DIRECTORY="integration-rci-${TEST_CONF}"
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
cp "${TEST_DIR}/isodata" isodata || exit 1
cp "${TEST_DIR}/${TEST_CONF}.c" "${TEST_CONF}.c" || exit 1

# Create the hydrogenic orbitals file (for W)
${GRASP_EXPORTHYDROGENIC} && mv "hydrogenic.w" "${TEST_CONF}.w" || {
	>&2 echo "ERROR: Failed to generate hydrogenic orbitals."
	exit 1
}

# Strip comments from rci.input
sed 's/[ \t]*#.*//' ${TEST_DIR}/rci.input > rci.input

# Using a temporary directory in order to not run into the bug where GRASP limits
# the length of the MPI_TMP variable.
export MPI_TMP=`mktemp -d`
>&2 echo "INFO: MPI_TMP=${MPI_TMP}"
>&2 echo "INFO: Running ${RCIQED}"
mpirun -n ${GRASP_MPI_NPROCS} --output-filename stdout ${RCIQED} < rci.input || {
	>&2 echo "FATAL ERROR: rci-qed failed with $?"
	>&2 echo "INFO: Cleanup:" && rm -Rv "${MPI_TMP}"
	exit 1
}
>&2 echo "INFO: Cleanup:" && rm -Rv "${MPI_TMP}"
>&2 echo "${RCIQED} ran successfully!"

# Verify the ouput
>&2 echo "INFO: Diffing ${TEST_CONF}.csum"
diff ${TEST_DIR}/reference.${TEST_CONF}.csum ${TEST_CONF}.csum || {
	>&2 echo "ERROR: There are differences in ${TEST_CONF}.csum"
	exit 1
}

# TODO: Also diff stdout. But it has (1) timing information, which is not
# deterministic, and (2) the MPI temporary directory is non-deterministic at
# the moment.
