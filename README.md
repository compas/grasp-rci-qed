# `rci-qed` program for GRASP

[![Doxygen Documentation](https://img.shields.io/badge/documentation-master-blue.svg)](http://mortenpi.eu/grasp-rci-qed/)
[![Build Status](https://travis-ci.com/mortenpi/grasp-rci-qed.svg?token=J2TJDmxGV6c9f8C3LXps&branch=master)](https://travis-ci.com/mortenpi/grasp-rci-qed)

**An updated RCI program with new QED for GRASP.**

This is a forked version of `rci90_mpi` from the 2018 version of GRASP that
includes an updated treatment of QED effects.

It needs to be linked against the `libmod`, `lib9290`, `libdvd90`, `librang90`
and `libmpiu90` libraries from GRASP. It also need the dependencies of these
libraries -- BLAS, LAPACK and MPI.

## Installation

For compiling it together with the published [2018 CPC version of
GRASP][compas-grasp2018], you should set the `GRASP` environment variable to
point to the GRASP root directory and call `./configure.sh` without any
arguments:

```
export GRASP=/path/to/grasp
./configure.sh
```

To integrate with the updated CMake-based build system of [the development
fork][mortenpi-grasp], source GRASP's `envset.sh` and pass the `--grasp-cmake`
argument to `./configure.sh`:

```
source /path/to/grasp/envset.sh
./configure.sh --grasp-cmake
```


[compas-grasp2018]: https://github.com/compas/grasp2018
[mortenpi-grasp]: https://github.com/mortenpi/grasp
