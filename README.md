# `rci-qed` program for GRASP

[![Doxygen Documentation][doc-img]][doc-url]
[![Build Status][travis-img]][travis-url]

**An updated RCI program with new QED for GRASP.**

This is a forked version of `rci90_mpi` from the 2018 version of GRASP that
includes an updated treatment of QED effects.

It needs to be linked against the `libmod`, `lib9290`, `libdvd90`, `librang90`
and `libmpiu90` libraries from GRASP. It also need the dependencies of these
libraries -- BLAS, LAPACK and MPI.

## Installation

For compiling it together with the published [2018 CPC version of GRASP][compas-grasp-2018],
you should set the `GRASP` environment variable to point to the GRASP root directory and
call `./configure.sh` without any arguments:

```
export GRASP=/path/to/grasp
./configure.sh
```

To integrate with the updated CMake-based build system available on [the GRASP master
branch][compas-grasp], source GRASP's `envset.sh` and pass the `--grasp-installed-modules`
argument to `./configure.sh`:

```
source /path/to/grasp/envset.sh
./configure.sh --grasp-installed-modules
```


[compas-grasp]: https://github.com/compas/grasp
[compas-grasp-2018]: https://github.com/compas/grasp/releases/tag/2018-12-03
[doc-img]: https://img.shields.io/badge/documentation-master-blue.svg
[doc-url]: http://compas.github.io/grasp-rci-qed/
[travis-img]: https://travis-ci.com/compas/grasp-rci-qed.svg?token=J2TJDmxGV6c9f8C3LXps&branch=master
[travis-url]: https://travis-ci.com/compas/grasp-rci-qed
