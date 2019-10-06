#!/bin/bash
if [ -z ${GRASP+x} ]; then >&2 echo "ERROR: \$GRASP variable is unset."; exit 1; fi

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TMPDIR="${DIR}/data/tmp"
DATADIR="${DIR}/data/belike"

source "${DIR}/../lib/testing.sh" || { echo "ERROR: Unable to load lib/testing.sh"; exit 1; }

# Locate the necessary GRASP binaries
RNUCLEUS=`find-grasp-binary rnucleus` || exit
RCSFGENERATE=`find-grasp-binary rcsfgenerate` || exit
RWFNESTIMATE=`find-grasp-binary rwfnestimate` || exit
RCIQED=`find-grasp-binary rci-qed` || exit
RCIQEDPT=`find-grasp-binary rci-qed.pt` || exit
echo "INFO: GRASP binary paths"
echo "  RNUCLEUS     = ${RNUCLEUS}"
echo "  RCSFGENERATE = ${RCSFGENERATE}"
echo "  RWFNESTIMATE = ${RWFNESTIMATE}"
echo "  RCIQED       = ${RCIQED}"
echo "  RCIQEDPT     = ${RCIQEDPT}"

if [ -d "${TMPDIR}" ]; then
	rm -rv "${TMPDIR}"
fi
mkdir -p "${TMPDIR}"
cd "${TMPDIR}" || exit 1

# Create the nuclear data file for a point-like xenon nucleus:
${RNUCLEUS} <<-EOF || { >&2 echo "ERROR: rnucleus failed with error code $?"; exit 1; }
	54
	0
	0
	0
	0
	0
EOF

# Generate CSLs for Be-like
${RCSFGENERATE} <<-EOF || { >&2 echo "ERROR: rcsfgenerate failed with error code $?"; exit 1; }
	*
	0
	1s(2,*)2s(2,*)

	3s,3p,3d
	0, 4
	4
	n
EOF
mv rcsf.out rcsf.inp || exit 1

${RWFNESTIMATE} <<-EOF || { >&2 echo "ERROR: rwfnestimate failed with error code $?"; exit 1; }
	y
	3
	*
EOF

# Set up hydrogen for CI
mv rwfn.inp belike.w || exit 1
mv rcsf.inp belike.c || exit 1

${RCIQED} <<-EOF || { >&2 echo "ERROR: rci-qed failed with error code $?"; exit 1; }
	y
	belike
	n
	n
	n
	n
	n
	n
	1-10
	1-9
	1-11
EOF

if ! [ -f "belike.cm" ]; then >&2 echo "ERROR: belike.cm is missing"; exit 1; fi

echo "Copying data over to: ${DATADIR}"
if [ -d "${DATADIR}" ]; then rm -rv "${DATADIR}"; fi
mkdir -p "${DATADIR}" || exit 1
cp isodata "${DATADIR}/isodata" || exit 1
cp belike.w "${DATADIR}/belike.w" || exit 1
cp belike.c "${DATADIR}/belike.c" || exit 1
cp belike.cm "${DATADIR}/belike.cm" || exit 1
