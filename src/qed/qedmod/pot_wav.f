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
! Originally from: f/potuse/pot_wav.f
!
c       =================================================
        subroutine SHAB_se_pot_wav(n,kappa,r,pp,qq,cp,cq)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the functions cp(i),cq(i),
c       which are the result of model potentail acting on
c       the wave functions p(i) and q(i) respectively.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_const/cl,z
        common /SHAB_grd/h,r1,r2,al,bt,ii
        common /SHAB_proj/ns_proj,num_kappa,nn_proj(maxns),
     &  ll_proj(maxns),jj_proj(maxns),kk_proj(maxns)
        common /SHAB_pots/v(maxii),unuc(maxii),uehl(maxii),v_wk(maxii),
     &  vloc(maxii),vsemi_loc(maxii,6),wsemi_loc(maxii,6),ww(maxii)
        common /SHAB_se/kk_se(6),ll_se(6),nn_se(6),se_loc(6),
     &  num_loc(-3:3),nshell(-5:5,10)
        common /SHAB_nloc/Dab(maxns,maxns),Gab(maxns,maxns),
     &  sab(maxns,maxns),vab(maxns,maxns),qed_matrix(maxns,maxns),
     &  de_qed0(maxns)
        common /SHAB_upuq/up(maxii,maxns),uq(maxii,maxns)
        real*8 r(maxii)
        real*8 pp(maxii),qq(maxii)
        real*8 cp(maxii),cq(maxii)
        real*8 dproj(maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        gam=dsqrt(kappa**2-(z/cl)**2)
        pp(ii+3)=ii
        qq(ii+3)=ii
        pp(ii+4)=gam
        qq(ii+4)=gam
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Local potential
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        il=0
        do i=1,num_kappa
          if (kappa.eq.kk_se(i)) then
            il=i
          endif
        enddo
        if (il.gt.0) then
          do i=1,ii
            cp(i)=vsemi_loc(i,il)*pp(i)/r(i)
            cq(i)=vsemi_loc(i,il)*qq(i)/r(i)
          enddo
        else
          do i=1,maxii
            cp(i)=0.d0
            cq(i)=0.d0
          enddo
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       non-local
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do nj=1, ns_proj
          dproj(nj)=0.d0
          if (kappa.eq.kk_proj(nj)) then
            dproj(nj)=SHAB_tint(-1,pp,qq,up(1,nj),uq(1,nj),r,v)
          endif
        enddo
c
        do nj=1, ns_proj
          dnonloc=0.d0
          do nk=1, ns_proj
            dnonloc=dnonloc+Dab(nj,nk)*dproj(nk)
          enddo
          do i=1,ii
            cp(i)=cp(i)+dnonloc*up(i,nj)/r(i)
            cq(i)=cq(i)+dnonloc*uq(i,nj)/r(i)
          enddo
        enddo
c
        cp(ii+3)=ii
        cq(ii+3)=ii
        cp(ii+4)=vsemi_loc(ii+4,il)+gam
        cq(ii+4)=vsemi_loc(ii+4,il)+gam
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
