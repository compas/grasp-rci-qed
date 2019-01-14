! This code is from QEDMOD with minor modifications.
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
c       =================================================
        subroutine SHAB_grid_args(r1_arg, r2_arg)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initializes the semi-logarithmic
c       radial grid r(i), which is used for the calculation
c       of the wave functions and the numerical integrations.
c
c       The nonuniform grid r is obtained from the uniform grid rho
c       by using the change of variables
c       rho(r)=al*r+bt*ln(r),
c       where al and bt are the grid parameters.
c
c       rho(r) is the uninform grid,
c       rho_i=rho_1+h*(i-1), where h is the grid spacing.
c
c       ii   .... number of grid points (ii=maxii-20), where
c                 maxii is given in the include file
c                 'qedmod.inc'
c       r(i) .... non-unform radial grid (i=1,ii)
c       r1   .... the first grid point closest to the origin
c       r2   .... the last grid point (the size of the box), defined
c                 in the atom_data subroutine r2=r(ii).
c       h    .... Semi-logarithmic grid step
c       z    .... Nuclear charge
c       rnucl ... Radius of the nuclear sphere
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_const/cl,z
        common /SHAB_grd/h,r1,r2,al,bt,ii
        common /SHAB_ncls/knucl,inucl,rnucl
        common /SHAB_pots/v(maxii),unuc(maxii),uehl(maxii),v_wk(maxii),
     &  vloc(maxii),vsemi_loc(maxii,6),wsemi_loc(maxii,6),ww(maxii)
        common /SHAB_ry/y(maxii),r(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r1=r1_arg
        r2=r2_arg
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        rmax=dabs(r2)
        imax=maxii-20
        imax=((imax+1)/2)*2-1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Semi-logarithmic grid: ro=al*r+bt*ln(r)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        ii=imax
        bt=1.d0
        h=bt*dlog(r2/r1)/((ii-1)*0.67)
        al=(h*(imax-1)-bt*dlog(r2/r1))/(r2-r1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (al.ge.0.d0) goto 210
        ii=(bt*dlog(r2/r1))/h+1
        write( *,99) imax,ii,h,al,bt,r1,r2
99      format(/'  ii > Imax'/'  Imax=',i6/'  ii  =',i6/
     1  '  h   =',f9.4,/'  al  =',f9.3/,'  bt  =',f9.3/
     2  '  r1  =',e10.3,'  r2  =',f9.3)
        call SHAB_exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Correction for the Sphere nuclear model
c       r(inucl)=rnucl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
210     inucl=0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (knucl.eq.2) then
          ro1=al*r1+bt*dlog(r1)
          ron=al*rnucl+bt*dlog(rnucl)
          d=ron-ro1
          inucl=(d/h+0.0001d0)+1
          if (inucl.lt.11) then
            inucl=11
            dr=d-(inucl-1)*h
            ro1=ro1+dr
            r0=dexp(ro1)/bt
            r1=SHAB_rr(ro1,r0)
          endif
          dr=(rnucl-r1)/(r2-r1)
          h=(bt*dlog(rnucl/r1)-bt*dlog(r2/r1)*dr)/
     &    ((inucl-1)-(ii-1)*dr)
          al=((inucl-1)*h-bt*dlog(rnucl/r1))/(rnucl-r1)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          r(i)=0.d0
          v(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call SHAB_tbr(ii,r1,h,al,bt,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        write( *,5) ii,r1,r2
        write(11,5) ii,r1,r2
5       format (/2x,'Semi-logarithmic radial grid:',
     1  /2x,'npoints    =',i6,/2x,'rmin    =',e13.4,
     2  /2x,'rmax    =',e13.4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
c       =================================================
        subroutine SHAB_tbr(ii,r1,h,al,bt,r,v)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Initialization of the Semi-logarithmic radial grid
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension r(*),v(*)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r(1)=r1
        d=r1
        t=al*d+bt
        v(1)=d/t
        p=al*d+bt*dlog(d)
        do 10 i=2,ii
        p=p+h
200     t=al*d+bt*dlog(d)
        t=(p-t)/(al*d+bt)
        d=d*(1.d0+t)
        if (dabs(t).gt.0.5d-11) goto 200
        t=al*d+bt
        v(i)=d/t
        r(i)=d
10      continue
        return
        end
c       =================================================
        function SHAB_rr(ro,r0)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine solve the equation ro=al*r+bt*ln(r)
c       and define the value of r for the given ro.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_grd/h,r1,r2,al,bt,ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        r=r0
200     t=al*r+bt*dlog(r)
        t=(ro-t)/(al*r+bt)
        if (t.le.-1.d0) t=-0.5d0
        r=r*(1.d0+t)
        if (dabs(t).gt.0.1d-8) goto 200
        SHAB_rr=r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
