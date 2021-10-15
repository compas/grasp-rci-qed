# Development notes

Replacing the old `rci90` and `rci90_mpi` programs with this version has the
following implications on the GRASP libraries:

* `where_C` in `libmod` is deprecated and should be removed
* Dependency on `ncdist_C` from `libmod` has been removed.
