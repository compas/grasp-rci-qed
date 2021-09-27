#!/usr/bin/env bash
#
# ./configure.sh [--debug] [--grasp-cmake]
#
# Configures a CMake-based out-of-tree build of grasp-rci-qed.
#
# By default, the build will be located under the build/ directory and will be 'Release'
# build. However, if `--debug` is passed, the build will instead be located under
# build-debug/ and have debug flags enabled. It is possible to have both build
# configurations be present in parallel.
#
# The `--grasp-cmake` option should be passed when GRASP itself was built using CMake. In
# that case, the `.mod` files get installed under ${GRASP}/lib/${library}, rather than next
# to the Fortran files under the `src/` directory. If this is passed, it properly passes the
# GRASP_INSTALLED_MODULES option to CMake, which makes sure that the build can find the
# `.mod` files that have been installed into `lib/` in the GRASP build.
#
# The script requires the $GRASP environment variable to be set, which must point to the
# root of a (make-installed) GRASP build. The CMake build for grasp-rci-qed will use the
# `.a` object archive and `.mod` files from the GRASP build.
#
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

# Expand GRASP into an absolute path, if not already:
export GRASP=$(realpath "$GRASP")

unset grasp_cmake_build
build_directory="build" # we default to build/
cmake_args=""
for arg in $@; do
	if [ "$arg" == "--debug" ]; then
		echo "Creating a DEBUG build"
		build_directory="build-debug"
		cmake_args="-DCMAKE_BUILD_TYPE=Debug${cmake_args:+ $cmake_args}"
	elif [ "$arg" == "--grasp-cmake" ]; then
		cmake_args="-DGRASP_INSTALLED_MODULES=true${cmake_args:+ $cmake_args}"
	fi
done

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

Note that you probably want to also enable parallel builds by passing -j to
the make command:

    make -jN

where N is the number of cores you have available.
EOF
