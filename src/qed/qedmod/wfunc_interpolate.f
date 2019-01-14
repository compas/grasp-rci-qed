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
! Originally from: f/potgen/potuse2.f
!
c=======================================================================
        subroutine SHAB_wfunc_interpolate(KA,NnA,x,GA,FA)
c
c 5-point interpolation of the stored wave function
c
        implicit real*8 (a-h,o-z)
        parameter (numRadPo = 10000)
        real*8 Fg(0:10),Ff(0:10)

        parameter (numwfMax = 10)
        common /SHAB_DF/xr(numRadPo,numwfMax),Gc(numRadPo,numwfMax),
     &  Fc(numRadPo,numwfMax),npoints(numwfMax),nc,KC(numwfMax),
     &  NnC(numwfMax)

        save left
        data left/0/
c--
        inum = 0
        do j = 1,nc
           if ((KA.eq.KC(j)).and.(NnA.eq.NnC(j))) inum = j
        enddo
        if (inum.eq.0) then
           print *,"kappa = ",KA," n = ",NnA
           stop 'wf_DF_int: 1'
        end if
c-
        GA = 0.d0
        FA = 0.d0
        nmax = npoints(inum)

        if (x.gt.xr(nmax,inum)) return
c-
        if ((left.ge.nmax).or.(left.lt.1)) left=1
        if (x.lt.xr(left,inum)) left=1

        do while ((left.lt.nmax).and.(x.gt.xr(left+1,inum)))
           left = left+1
        enddo

        imin = max(1,   left-2)
        imax = min(nmax,left+2)
        if (imin.eq.1) imax = 3
        if (imax.eq.nmax) imin = nmax-4

        do j = imin,imax
           Fg(j-imin) = Gc(j,inum)
           Ff(j-imin) = Fc(j,inum)
        enddo

        do j = imin+1,imax
           do i = j,imax
              Fg(i-imin) = (Fg(j-imin-1)-Fg(i-imin))
     >                                /(xr(j-1,inum)-xr(i,inum))
              Ff(i-imin) = (Ff(j-imin-1)-Ff(i-imin))
     >                                /(xr(j-1,inum)-xr(i,inum))
           enddo
        enddo

        GA = Fg(imax-imin)
        FA = Ff(imax-imin)
        do i = imax-1,imin,-1
           GA = Fg(i-imin)+ (x-xr(i,inum))* GA
           FA = Ff(i-imin)+ (x-xr(i,inum))* FA
        enddo
c--
        return
        end
