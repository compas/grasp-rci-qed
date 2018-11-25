# Testing-related routines

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

function check_files {
	for filename in "$@"; do
		if ! [ -f "${filename}" ]; then
			echo "ERROR: ${filename} missing."
			exit 2
		fi
	done
}

function check_directories {
	for dirname in "$@"; do
		if ! [ -d "${dirname}" ]; then
			echo "ERROR: ${dirname} does not exist."
			exit 2
		fi
	done
}
