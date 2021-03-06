#!/bin/bash
export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "${SCRIPT_DIR}/../lib/testing.sh" || { echo "ERROR: Unable to load lib/testing.sh"; exit 1; }

if [ -z ${GRASP+x} ]; then >&2 echo "ERROR: \$GRASP variable is unset."; exit 1; fi
if [ -z ${GRASP_BUILD_BINDIR+x} ]; then >&2 echo "ERROR: \$GRASP_BUILD_BINDIR variable is unset."; exit 1; fi
if [ -z ${GRASP_EXPORTHYDROGENIC+x} ]; then >&2 echo "ERROR: \$GRASP_EXPORTHYDROGENIC variable is unset."; exit 1; fi
if ! [ -f ${GRASP_EXPORTHYDROGENIC} ]; then >&2 echo "ERROR: \$GRASP_EXPORTHYDROGENIC not pointing to a file."; exit 1; fi

# We'll need to locate the rci-qed and rci-qed.pt binaries
RCIQED=`find-grasp-binary rci-qed` || exit
echo "INFO: RCIQED=${RCIQED}"
RCIQEDPT=`find-grasp-binary rci-qed.pt` || exit
echo "INFO: RCIQEDPT=${RCIQEDPT}"
RCIQEDORBITALS=`find-grasp-binary rci-qed.orbitals` || exit
echo "INFO: RCIQEDORBITALS=${RCIQEDORBITALS}"

# Get the configuration -- the first argument, and check that input files exist
if [ "$1" == "--update" ]; then
	>&2 echo "INFO: Update mode -- updating reference values."
	UPDATE_MODE=true
	shift
fi

TEST_CONF=$1
TEST_DIR="${SCRIPT_DIR}/$1"
check_directories "${TEST_DIR}"
check_files "${TEST_DIR}/isodata" "${TEST_DIR}/rci.input" "${TEST_DIR}/test.c"

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
cp "${TEST_DIR}/test.c" "test.c" || exit 1

# Create the hydrogenic orbitals file (for W)
${GRASP_EXPORTHYDROGENIC} && mv "hydrogenic.w" "test.w" || {
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

# Run rci-qed.pt
>&2 echo "INFO: Running ${RCIQEDPT}"
${RCIQEDPT} test 2>&1 | tee rci-qed.pt.stdout
if ! [ "${PIPESTATUS[0]}" == "0" ]; then
	>&2 echo "FATAL ERROR: rci-qed.pt failed with ${PIPESTATUS[0]}"
	exit 1
fi
>&2 echo "${RCIQEDPT} ran successfully!"

# Run rci-qed.orbitals
>&2 echo "INFO: Running ${RCIQEDORBITALS}"
${RCIQEDORBITALS} test 2>&1 | tee rci-qed.orbitals.stdout
if ! [ "${PIPESTATUS[0]}" == "0" ]; then
	>&2 echo "FATAL ERROR: rci-qed.orbitals failed with ${PIPESTATUS[0]}"
	exit 1
fi
>&2 echo "${RCIQEDORBITALS} ran successfully!"

# Verify the ouput
function verify_file {
	filename=$1
	>&2 echo "INFO: Diffing ${filename}"
	diff "${TEST_DIR}/${filename}.reference" "${filename}" || {
		>&2 echo "ERROR: There are differences in ${filename}"
		if [ "$UPDATE_MODE" == "true" ]; then
			>&2 echo "INFO: updating ${TEST_DIR}/${filename}.reference"
			mv "${filename}" "${TEST_DIR}/${filename}.reference"
		else
			exit 1
		fi
	}
}

verify_file "test.csum"
verify_file "test.settings.toml"
verify_file "rci-qed.pt.stdout"
verify_file "rci-qed.orbitals.stdout"


# TODO: Also diff stdout. But it has (1) timing information, which is not
# deterministic, and (2) the MPI temporary directory is non-deterministic at
# the moment.
