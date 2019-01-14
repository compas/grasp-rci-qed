# QED subroutines

Contains routines for calculating QED self-energies using different approaches.
The routines are compiled into a static `libqed.a` library. It directly depends
on routines and common blocks from `lib92`.

The `legacy/` subdirectory contains code that has been moved over from the `rci`
application with no or minimal changes. Those routines are necessary for the
functioning of the new routines, but they should be refactored sooner or later.


## Modified `QEDMOD` routines

Fortran routines for QED self-enery and vacuum polarization (Uehling, Wichmann-Kroll)
calculations using the model operator approach live under the `shabaev/` directory.
Code originally distributed in CPC:

> V.M. Shabaev, I.I. Tupitsyn, V.A. Yerokhin, QEDMOD: Fortran program for calculating
> the model Lamb-shift operator, Computer Physics Communications, Volume 189, 2015,
> Pages 175-181, ISSN 0010-4655, http://dx.doi.org/10.1016/j.cpc.2014.12.002.

### Modifications

The code redistributed here contains a subset of the routines and is slightly modified.

* All routine and function names were prepended with `SHAB_` so that they would not
  clash with existing routines when linking against other software.

* Common blocks were clumped into (somewhat) logical groups and renamed using the
  `SHAB_` prefix in order to prevent clashes with the Grasp2k code.

* Various routines were split into separate files:

  - `exit` into `exit.f`.
  - `init_se` from `potgen/se_pot.f` because it doesn't depend on any of the
    other routines in that file and the other routines back-depend on
    `potgen_main.f`.
  - `read_func` from `potgen/potgen_main.f` to `read_func.f`, since the routines
    in `se_pot.f` depend on it.
  - break out the piece in `potgen_main` that populates the `/ff/` common block with
    hydrogenic wavefunctions. This is now a separate routine
    (`SHAB_populate_hydrogenics`) in `populate_hydrogenics.f90`
  - move `write_func` from `potgen/potgen_main.f` into `write_func.f`, because it
    is used by `SHAB_populate_hydrogenics`.
  - move `wfunc_interpolate` from `potuse/potuse2.f` to `wfunc_interpolate.f`,
    so the that interpolation routine could be used by external code.

* Other routines were removed from their corresponding files since they are not
  relevant for the implementation in Grasp2k:

  - `write_pot` in `potgen/se_pot.f`, since only `potgen_main.f` calls that routine.

* Add `set_fermi_params` argument to the `nucl` routine. If set to `.false.`,
  then the `tt`, `aa` and `cc` variables are not set by `nucl` and are assumed
  to be set beforehand by the user. This is necessary to pass those values from
  Grasp2k to this code.


## Credits

* The original cleanup of some of the existing QED routines and implementation of
  them as a separete library was done by Lukáš Félix Pašteka and Morten Piibeleht
  in 2017.
* The Flambaum and Pyykkö implementation was based off an implementation of these
  routines in Grasp92 by Christian Thierfelder, with modifications by Peter
  Schwerdtfeger and Lukáš Félix Pašteka.
* We're also redistributing a subset of modified routines from
  [QEDMOD by Shabaev et al.](http://dx.doi.org/10.1016/j.cpc.2014.12.002), which
  can be found under the `shabaev/` directory.
* Various routines were stripped out of the `rci` program and moved here with
  minor or no changes and can be found under `legacy/`.
