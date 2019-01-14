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
        subroutine SHAB_local_se_pot
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the local part Vloc(r)
c       of the Self-Energy (SE) QED potential.
c       Local potential multiplied on r is store in the
c       array 'vsemi_loc(i,kappa)' for each value of the
c       relativistic quantum number kappa, where 'i' is
c       a number of grid point.
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
        common /SHAB_ry/y(maxii),r(maxii)
        common /SHAB_pq/p(maxii),q(maxii)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do il=1,num_kappa
          do i=1,maxii
            vsemi_loc(i,il)=0.d0
            wsemi_loc(i,il)=0.d0
          enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 il=1,num_kappa
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_se(il)
        l=ll_se(il)
        n=nn_se(il)
        ni=nshell(kappa,n)
        if (ni.eq.0) goto 10
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          ri=r(i)
          vloc(i)=dexp(-ri*cl)*ri
          ww(i)=  dexp(-ri*cl)*ri
        enddo
        do i=ii,1,-1
          if (dabs(vloc(i)).gt.1.d-40) goto 200
          imax=i-1
          vloc(i)=0.d0
        enddo
c
 200    vloc(ii+3)=imax
        vloc(ii+4)=1.d0
        ww(ii+3)=ii
        ww(ii+4)=1.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call SHAB_read_func(ni,p,q,2)
        de=SHAB_sint(vloc,p,q,p,q,r,v)
        vloc0=se_loc(il)/de
        de=SHAB_sint(ww,p,q,p,q,r,v)
        w0=1.d0/de
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,ii
          vloc(i)=vloc0*vloc(i)
          ww(i)=w0*ww(i)
        enddo
        do i=1,maxii
          vsemi_loc(i,il)=vloc(i)
          wsemi_loc(i,il)=ww(i)
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine SHAB_nonlocal_se_pot
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes the nonlocal (separable)
c       part of the Self-Energy (SE) potential. This
c       potential has the form:
c       V_nonloc = \sum_{a,b} |a> D_ab <b|, where a,b
c       are projected wave functions and
c       Dab=Sab^{-1} * QED_matr * Sab^{-1}
c       The projected wave functions and matrix D_ab are
c       stored in the external file
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
        common /SHAB_ry/y(maxii),r(maxii)
        common /SHAB_pq/p(maxii),q(maxii)
        common /SHAB_cpcq/cp(maxii),cq(maxii)
        common /SHAB_upuq/up(maxii,maxns),uq(maxii,maxns)
        real*8 a(maxii),b(maxii)
        real*8 de_qed_loc0(maxns)
        real*8 diag(maxns),diag_qed(maxns),diag_inv(maxns)
        real*8 p1(maxii),q1(maxii)
        real*8 p2(maxii),q2(maxii)
        real*8 p3(maxii),q3(maxii)
        integer iwrk1(maxns),iwrk2(maxns)
        parameter (np=5)
        character*1 let(11)
        data
     1  let /'s','p','d','f','g','h','i','k','l','m','n'/
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter (pi=3.1415926535897932385d0)
        alpha=1.d0/cl
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Functions up, uq
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
          kappa=kk_proj(ni)
          n=nn_proj(ni)
          l=ll_proj(ni)
          call SHAB_read_func(ni,p,q,2)
          il=num_loc(kappa)
          do i=1,maxii
            vloc(i)=vsemi_loc(i,il)
            ww(i)=wsemi_loc(i,il)
          enddo
          do i=1,maxii
            a(i)=p(i)
            b(i)=q(i)
          enddo
          do i=1,ii
            a(i)=ww(i)*p(i)
            b(i)=ww(i)*q(i)
          enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          k=nn_proj(ni)-ll_proj(ni)
            do i=1,ii
              ppi=p(i)
              qqi=q(i)
              dw=dexp(-r(i)*z*2)*r(i)
              factor_p=(1.d0-(-1.d0)**k)/2.d0
              factor_q=(1.d0+(-1.d0)**k)/2.d0
              a(i)= factor_p*dw*ppi
              b(i)= factor_q*dw*qqi
            enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
          gam=p(ii+4)+ww(ii+4)
          a(ii+4)=gam
          b(ii+4)=gam
          a(ii+3)=ii
          b(ii+3)=ii
c
          dn=SHAB_tint(-1,p,q,a,b,r,v)
          do i=1,ii
            up(i,ni)=a(i)/dn
            uq(i,ni)=b(i)/dn
          enddo
          imax=ii
          do i=ii,1,-1
            if (dabs(up(i,ni))+dabs(uq(i,ni)).gt.1.d-40) goto 200
            imax=i-1
            up(i,ni)=0.d0
            uq(i,ni)=0.d0
          enddo
 200      up(ii+3,ni)=imax
          uq(ii+3,ni)=imax
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Sab and Vab
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kappa=kk_proj(ni)
        il=num_loc(kappa)
        do i=1,maxii
          vloc(i)=vsemi_loc(i,il)
          ww(i)=wsemi_loc(i,il)
        enddo
        call SHAB_read_func(ni,p,q,2)
        do nj=1,ns_proj
          sab(ni,nj)=0.d0
          vab(ni,nj)=0.d0
          if (kk_proj(ni).eq.kk_proj(nj)) then
            call SHAB_read_func(nj,cp,cq,2)
            de=SHAB_sint(vloc,p,q,cp,cq,r,v)
            vab(ni,nj)=de
            do i=1,maxii
              a(i)=up(i,nj)
              b(i)=uq(i,nj)
            enddo
            sab(ni,nj)=SHAB_tint(-1,p,q,a,b,r,v)
          endif
        enddo
        de_qed_loc0(ni)=vab(ni,ni)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Matrix Inversion
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          Gab(ni,nj)=Sab(ni,nj)
        enddo
        enddo
        call SHAB_dminv(Gab,ns_proj,maxns,det,iwrk1,iwrk2)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Inversion test
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        del=0.d0
        do ni=1,ns_proj
        do nj=1,ns_proj
          aij=0.d0
          do nk=1,ns_proj
            aij=aij+Gab(ni,nk)*Sab(nk,nj)
          enddo
          if (ni.eq.nj) aij=aij-1.d0
          del=del+dabs(aij)/ns_proj
        enddo
        enddo
        if (abs(del).lt.1.d-10) then
           write( *,'(/2x,a,e12.2,a)') 'Inversion test:',del,"  OK"
           write(11,'(/2x,a,e12.2,a)') 'Inversion test:',del,"  OK"
        else
           write( *,'(/2x,a,e12.2,a)') 'Inversion test:',del
           write(11,'(/2x,a,e12.2,a)') 'Inversion test:',del
        end if
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call SHAB_dab_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1, ns_proj
          delt=0.d0
          do nj=1, ns_proj
          do nk=1, ns_proj
            delt=delt+Sab(ni,nj)*Dab(nj,nk)*Sab(ni,nk)
          enddo
          enddo
          de_qed0(ni)=Vab(ni,ni)+delt
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 1000   return
        end
c       =================================================
        subroutine SHAB_dab_matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       The subroutine computes Dab matrix
c       Dab=Sab^{-1} * QED_matr * Sab^{-1}
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        include 'qedmod.inc'
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_proj/ns_proj,num_kappa,nn_proj(maxns),
     &  ll_proj(maxns),jj_proj(maxns),kk_proj(maxns)
        common /SHAB_nloc/Dab(maxns,maxns),Gab(maxns,maxns),
     &  sab(maxns,maxns),vab(maxns,maxns),qed_matrix(maxns,maxns),
     &  de_qed0(maxns)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       D-matrix
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do ni=1,ns_proj
        do nj=1,ns_proj
          dab(ni,nj)=0.d0
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 10 ni=1,ns_proj
        do 20 nj=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        Dab(ni,nj)=0.d0
        if (kk_proj(nj).ne.kk_proj(ni)) goto 20
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 30 nk=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nk).ne.kk_proj(ni)) goto 30
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do 40 nl=1,ns_proj
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (kk_proj(nl).ne.kk_proj(ni)) goto 40
        qkl=qed_matrix(nk,nl)
        qlk=qed_matrix(nl,nk)
        dkl=0.5d0*(qkl+qlk)-vab(nk,nl)
        Dab(ni,nj)=Dab(ni,nj)+Gab(ni,nk)*dkl*Gab(nj,nl)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
40      continue
30      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
20      continue
10      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine SHAB_dminv(a,n,nmax,d,l,m)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Invertion of arbitrary matrix:  C=A**(-1),
c       using Gauss-Gordan method.
c       This subroutine is a sligtly modified version of the
c       minv subroutine from the IBM Application Program:
c       Scientific Subroutine Package (SSP)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       l,m - work arrays (dimension - n)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        dimension a(nmax,nmax),l(nmax),m(nmax)
c       double precision a,d,biga,hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Searching maximal element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=1.d0
        do 80 k=1,n
        l(k)=k
        m(k)=k
        biga=a(k,k)
        do 20 j=k,n
        do 20 i=k,n
10      if (dabs(biga)-dabs(a(i,j))) 15,20,20
15      biga=a(i,j)
        l(k)=i
        m(k)=j
20      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the strings
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        j=l(k)
        if (j-k) 35,35,25
25      do 30 i=1,n
        hold=-a(k,i)
        a(k,i)=a(j,i)
30      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Changing the columns
c       - - - - - - - - - - - - - - - - - - - - - - - - -
35      i=m(k)
        if (i-k) 45,45,38
38      jp=n*(i-1)
        do 40 j=1,n
        hold=-a(j,k)
        a(j,k)=a(j,i)
40      a(j,i)=hold
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding column on the leader element.
c       The leader element is in 'biga'.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
45      if (biga) 48,46,48
46      d=0.d0
        return
48      biga=-biga
        do 55 i=1,n
        if (i-k) 50,55,50
50      a(i,k)=a(i,k)/biga
55      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        biga=-biga
        do 65 i=1,n
        hold=a(i,k)
        ij=i-n
        do 65 j=1,n
        ij=ij+n
        if (i-k) 60,65,60
60      if (j-k) 62,65,62
62      kj=ij-i+k
        a(i,j)=hold*a(k,j)+a(i,j)
65      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Deviding the string on the leader element.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        kj=k-n
        do 75 j=1,n
        kj=kj+n
        if (j-k) 70,75,70
70      a(k,j)=a(k,j)/biga
75      continue
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Composition of the leader elements.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        d=d*biga
        a(k,k)=1.d0/biga
80      continue
        k=n
100     k=k-1
        if (k) 150,150,105
105     i=l(k)
        if (i-k) 120,120,108
108     do 110 j=1,n
        hold=a(j,k)
        a(j,k)=-a(j,i)
110     a(j,i)=hold
120     j=m(k)
        if (j-k) 100,100,125
125     do 130 i=1,n
        hold=a(k,i)
        a(k,i)=-a(j,i)
130     a(j,i)=hold
        goto 100
c       - - - - - - - - - - - - - - - - - - - - - - - - -
150     return
        end
