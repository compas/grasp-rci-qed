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
! Originally from: f/potgen/potgen_main.f
!
module SHAB_populate_hydrogenics_common
    implicit real*8 (a-h,o-z)

    include 'qedmod.inc'

    COMMON /SHAB_grd/h,r1,r2,al,bt,ii
    COMMON /SHAB_proj/ns_proj,num_kappa,nn_proj(maxns),ll_proj(maxns),  &
                      jj_proj(maxns),kk_proj(maxns)
    COMMON /SHAB_num/vnuc(20),nmax,niter,nit,m1,m2,m3

    character*1, parameter :: let(11) = (/'s','p','d','f','g','h','i','k','l','m','n'/)
end module SHAB_populate_hydrogenics_common

!> Populates the `/ff/` common block with hydrogenic wavefunctions.
!!
!! `verbose` is a boolean argument, and the routine will print stuff if set to true.
!! If `log_fhandle` is greater than 0, then the routine also writes the same stuff
!! into the file opened on that file handle. If set to 0 then nothing gets written.
!!
!! Adapted from a piece of the `potgen_main` program in `potgen/potgen_main.f`.
!!
subroutine SHAB_populate_hydrogenics(verbose, log_fhandle)
    use SHAB_populate_hydrogenics_common
    implicit none

    logical :: verbose
    integer :: log_fhandle

    character(*), parameter :: fmt = "(i3,i4,a1,i2,'/2',f16.8,i9)"

    integer :: ni, n, l, kappa
    real*8 :: d, e
    real*8 :: p(maxii), q(maxii)

    do ni = 1, ns_proj
        n=nn_proj(ni)
        l=ll_proj(ni)+1
        kappa=kk_proj(ni)
        call SHAB_dirac(n,kappa,p,q)
        e = 0.5d0*p(ii+1)
        d = p(ii+5)**2+q(ii+5)**2

        if(verbose) then
            write( *,fmt) ni,n,let(l),jj_proj(ni),e,niter
        end if
        if(log_fhandle > 0) then
            write(log_fhandle, fmt) ni,n,let(l),jj_proj(ni),e,niter
        end if

        call SHAB_write_func(ni,p,q,2)
    end do
end subroutine SHAB_populate_hydrogenics
