# Development notes

Replacing the old `rci90` and `rci90_mpi` programs with this version has the
following implications on the GRASP libraries:

* `where_C` in `libmod` is deprecated and should be removed
* Dependency on `ncdist_C` from `libmod` has been removed, which used to provide the `ZDIST` global variable storing the nuclear charge distributions and, later, the vacuum polarization potential.
* Dependency on `vpilst_C` from `libmod` has been removed, which was used to dynamically store the VP matrix elements for lookup. Relatedly, `NVPItmp` variable from `hmat_C` is no longer used.
