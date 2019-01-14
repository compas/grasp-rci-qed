! This code is adapted from QEDMOD.
!
!   V.M. Shabaev, I.I. Tupitsyn, V.A. Yerokhin, QEDMOD: Fortran program for calculating
!   the model Lamb-shift operator, Computer Physics Communications, Volume 189, 2015,
!   Pages 175-181, ISSN 0010-4655, http://dx.doi.org/10.1016/j.cpc.2014.12.002.
!
! For a description of the changes made, see the README.
!
! Original code is distributed under the CPC non-profit use licence
! (http://cpc.cs.qub.ac.uk/licence/licence.html). Redistributed here as part of
! GRASP2K with permission from the authors.
!
! Originally from: f/potgen/grid.f
!
module SHAB_grid_common
    implicit real*8 (a-h,o-z)

    include 'qedmod.inc'

    COMMON /SHAB_const/cl,z
    COMMON /SHAB_grd/h,r1,r2,al,bt,ii
    COMMON /SHAB_ncls/knucl,inucl,rnucl
end module SHAB_grid_common

!>
subroutine SHAB_grid
    use SHAB_grid_common
    implicit none

    real*8 :: r1_local, r2_local

    r1_local=dexp(-6.d0)*(256.d0/maxii)**2/5.d0
    if (z.gt.0.99) then
      r1_local=r1_local/z
    else
      r1_local=r1_local/100.d0
    endif
    if (rnucl > 0 .and. r1_local > 0.10d0*rnucl) r1_local=0.10d0*rnucl

    r2_local = r2

    call SHAB_grid_args(r1_local, r2_local)
end subroutine SHAB_grid
