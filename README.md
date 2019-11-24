# `rci-qed` for GRASP

[![Doxygen Documentation][doc-img]][doc-url]
[![Build Status][travis-img]][travis-url]

**An updated RCI program with new QED for GRASP.**

> **Beware:** this version is still a work on progress, contains unfinished features and probably bugs. Do not rely on this for scientific calculations.

This is a forked version of `rci90_mpi` from the 2018 version of GRASP that
includes an updated treatment of QED effects.

It needs to be linked against the `libmod`, `lib9290`, `libdvd90`, `librang90`
and `libmpiu90` libraries from GRASP. It also need the dependencies of these
libraries -- BLAS, LAPACK and MPI.

## Installation

To download the latest version of the code, you can clone it with Git:

```
git clone https://github.com/compas/grasp-rci-qed.git
```

When configuring `rci-qed` for compilation, you first need to make sure that you have the `$GRASP` environment variable pointing to the root of the GRASP directory:

```
export GRASP=/path/to/grasp
```

The compilation and installation scripts assume that the static libraries
are stored under `$GRASP/lib` and that compiled binaries should be installed to `$GRASP/bin`.

Next, the `./configure.sh` script can help setting up the [CMake](https://cmake.org/) build of `rci-qed`, by creating an _out-of-tree_ build directory under `build/`:

```
./configure.sh
```

After this, the `rci-qed` build directory is set up under `build/`.

### Compilation

* To compile the code and binaries, change into the `build/` directory and call `make`:

  ```
  cd build/
  make
  ```

* To run the test suite, call `ctest` or `make test` in `build`.

* To install the binaries into `$GRASP/bin`, run `make install` in `build/`.


[compas-grasp2018]: https://github.com/compas/grasp2018
[mortenpi-grasp]: https://github.com/mortenpi/grasp
[doc-img]: https://img.shields.io/badge/documentation-master-blue.svg
[doc-url]: http://compas.github.io/grasp-rci-qed/
[travis-img]: https://travis-ci.com/compas/grasp-rci-qed.svg?token=J2TJDmxGV6c9f8C3LXps&branch=master
[travis-url]: https://travis-ci.com/compas/grasp-rci-qed
