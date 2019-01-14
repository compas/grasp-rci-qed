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
! Originally from: f/potgen/se_pot.f
!
c       =================================================
        subroutine SHAB_init_se
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine initializes some data needed to
c       construct self-energy QED potential. In partuclar,
c       it computes the diagonal and non-diagonal SE matrix
c       elements by means of subroutine FSE_dat and stores them
c       in the array qed_matrix.
c
c       The quantum numbers listed below define the quantum
c       numbers of the projected wave functions used by the
c       nonlocal part of SE model potential.
c
c       Parameters of the projected wave functions:
c       ns_proj     - number of atomic shells
c       maxns       - maximal number of atomic shells
c                     (see file 'qedmod.inc')
c       ni          - shell number (ni=1,,,ns)
c       nn_proj(ni) - principal quantum number
c       ll_proj(ni) - orbital quantum number
c       jj_proj(ni) - angular quantum number
c       kk_proj(ni) - relativistic quantum numbers
c       num_kappa   - total number of different relativistic
c                     angular quantum numbers
c
c       z    .... Nuclear charge
c       cl   ..... Speed of light
c       knucl  ..... Nuclear model (0 - point, 1 - Fermi model,
c                    2 - uniform sphere)
c
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_const/cl,z
        common /SHAB_ncls/knucl,inucl,rnucl
        common /SHAB_proj/ns_proj,num_kappa,nn_proj(maxns),
     &  ll_proj(maxns),jj_proj(maxns),kk_proj(maxns)
        common /SHAB_se/kk_se(6),ll_se(6),nn_se(6),se_loc(6),
     &  num_loc(-3:3),nshell(-5:5,10)
        common /SHAB_nloc/Dab(maxns,maxns),Gab(maxns,maxns),
     &  sab(maxns,maxns),vab(maxns,maxns),qed_matrix(maxns,maxns),
     &  de_qed0(maxns)
        integer nn_max(0:3)
        parameter (np=5)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Default value of the max principal quantum number
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        nma=3
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do l=0,3
          if (nma.gt.0) then
            nn_max(l)=nma
          else
            nn_max(l)=iabs(nma)+l
          endif
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kk_se(1)=-1
        kk_se(2)=1
        kk_se(3)=-2
        kk_se(4)=2
        kk_se(5)=-3
        kk_se(6)=3
        num_kappa=5
        ni=0
        do il=1,num_kappa
          kappa=kk_se(il)
          l=(iabs(2*kappa+1)-1)/2
          ll_se(il)=l
          nn_se(il)=l+1
          do n=l+1,nn_max(l)
            ni=ni+1
            kk_proj(ni)=kappa
            ll_proj(ni)=l
            jj_proj(ni)=2*iabs(kappa)-1
            nn_proj(ni)=n
          enddo
        enddo
        ns_proj=ni
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        num_loc(-3)=5
        num_loc(-2)=3
        num_loc(-1)=1
        num_loc( 0)=0
        num_loc( 1)=2
        num_loc( 2)=4
        num_loc( 3)=6
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do kappa=-5,5
        do n=1,10
          nshell(kappa,n)=0
        enddo
        enddo
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          nshell(kappa,n)=ni
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation of the qed_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call SHAB_se_pnt_stor()
        call SHAB_se_fn_stor()
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=ni,ns_proj
          qed_matrix(ni,nj)=0.d0
          kappa=kk_proj(ni)
          nn1=nn_proj(ni)
          nn2=nn_proj(nj)
          if (kk_proj(nj).eq.kappa) then
            iz=z+0.5d0
            call SHAB_FSE_dat(kappa,nn1,nn2,iz,FSE_pnt,FSE_ext)
            dn1=nn1
            dn2=nn2
            dn=dsqrt(dn1*dn2)
            coef=alpha/pi*(z*alpha)**4/dn**3*cl**2
            if (knucl.eq.0) then
              qed_matrix(ni,nj)=fse_pnt*coef
              qed_matrix(nj,ni)=fse_pnt*coef
            else
              qed_matrix(ni,nj)=fse_ext*coef
              qed_matrix(nj,ni)=fse_ext*coef
            endif
            do il=1,num_kappa
              n=nn_se(il)
              if (kappa.eq.kk_se(il).and.nn1.eq.n.and.nn2.eq.n) then
                se_loc(il)=qed_matrix(ni,nj)
              endif
            enddo
          endif
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
