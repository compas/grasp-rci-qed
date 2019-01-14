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
! Originally from: f/potgen/potgen_main.f
!
c       =================================================
        subroutine SHAB_read_func(ni,v1,v2,nrec)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       reading a vector from the common block
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(maxff=2*(2*maxns+1))
        common /SHAB_ff/ff(maxii,maxff)
        integer ni,nrec
        real*8 v1(maxii),v2(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nr1=2*ni-1
        nr2=nr1+1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,maxii
          v1(i)=ff(i,nr1)
        enddo
        if (nrec.eq.1) goto 1000
        do i=1,maxii
          v2(i)=ff(i,nr2)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
1000    return
        end
