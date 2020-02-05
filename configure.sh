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
	if [ "$arg" == "--debug" ]; then
		echo "Creating a DEBUG build"
		build_directory="build-debug"
		cmake_args="-DCMAKE_BUILD_TYPE=Debug${cmake_args:+ $cmake_args}"
	elif [ "$arg" == "--grasp-installed-modules" ]; then
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
