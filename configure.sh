#!/usr/bin/env bash
if [ -z ${GRASP+x} ]; then
	>&2 cat <<-EOF
		ERROR: \$GRASP environment variable unset.
		It should point to the GRASP root directory. Set it with:

		    export GRASP=/path/to/grasp

		Or by sourcing the GRASP envset.sh script.
	EOF
	exit 1
fi
export DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


unset grasp_cmake_build
build_directory="build" # we default to build/
cmake_args=""
for arg in $@; do
	if [ "$arg" == "--grasp-cmake" ]; then
		grasp_cmake_build="enabled"
	fi
done
for arg in $@; do
	if [ "$arg" == "--debug" ]; then
		if ! [ -z ${grasp_cmake_build+x} ]; then
			>&2 echo "ERROR: can't specify --debug together with --grasp-cmake"
			exit 1
		fi
		echo "Creating a DEBUG build"
		build_directory="build-debug"
		cmake_args="-DCMAKE_BUILD_TYPE=Debug"
	elif [ "$arg" == "--grasp-cmake" ]; then
		unset build_directory
	fi
done

if ! [ -z ${grasp_cmake_build+x} ]; then
	# If --grasp-cmake was specified, we'll try to hook into the standard GRASP build.
	if ! [ -f "${GRASP}/CMakeLists.txt" ]; then
		>&2 echo "ERROR: CMakeLists.txt missing in ${GRASP}"
		exit 1
	fi
	cmakelistsuser="${GRASP}/CMakeLists.user"
	echo "Adding grasp-rci-qed to the CMakeLists.user file of GRASP"
	echo "  appending to: ${cmakelistsuser}"
	cat >> ${cmakelistsuser} <<-EOF
		# Added automatically by the ./configure.sh script of grasp-rci-qed
		add_subdirectory("${DIR}/src/" "\${CMAKE_CURRENT_BINARY_DIR}/external/grasp-rci-qed")
	EOF
	cat <<-EOF
		You can now build rci-qed when building GRASP.
		You may need to recreate the build/ directory.
		When you call \`make install\`, rci-qed will also be installed next to
		the other GRASP binaries in \${GRASP}/bin/.
	EOF
	exit
else
	# If not, we run the default setup for CMake's out-of-tree builds
	#
	#     mkdir build/
	#     cd build/
	#     cmake ..
	#
	build_abspath="${DIR}/${build_directory}"
	echo "Creating: ${build_abspath}"
	if [ -e "${build_abspath}" ]; then
		>&2 echo "ERROR: Build directory already exists."
		exit 1
	fi
	mkdir "${build_abspath}" && cd "${build_abspath}" \
		&& cmake ${cmake_args} "${DIR}"

	# Note: we need to use spaces, not tabs, to indent in the heredoc.
	cat <<-EOF
	Build directory ${build_directory}/ created.
	To build rci-qed run you need to cd into ${build_directory}/ and run make:

	    cd ${build_directory}/
	    make

	Note that you probably want to also enable parallel builds by pass -j to make:

	    make -jN

	where N is the number of cores you have available.
	EOF
fi
