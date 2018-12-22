#!/bin/bash
# Uses `diff` to verify whether two files are identical or not.
#
#   ./filematch.sh <output-file> <reference-file>
#
output=$1
reference=$2

if ! [ -f "${output}" ]; then
	>&2 echo "ERROR: Output file missing: ${ouput}"
	exit 2
fi
if ! [ -f "${reference}" ]; then
	>&2 echo "ERROR: Reference file missing: ${reference}"
	exit 2
fi

# Verify the ouput
>&2 echo "INFO: Diffing ${output}"
diff "${reference}" "${output}" || {
	>&2 echo "ERROR: There are differences in ${output}"
	if [ "$UPDATE_MODE" == "true" ]; then
		>&2 echo "INFO: updating ${reference}"
		mv "${output}" "${reference}"
	fi
	exit 1
}
