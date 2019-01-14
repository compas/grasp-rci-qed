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
! Originally from: f/potgen/wk.f
!
c       =================================================
        subroutine SHAB_wk(maxii,r,v_wk)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       This subroutine computes Wichmann-Kroll point nuclear
c       vacuum-polarization potential by means of
C       the analytical approximation formulas from
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569.
c
c       Input parameters:
c       r(i)    .... Radial grid
c       maxii   .... the size of the r array
c
c       Output parameters:
c       wk(i)   .... Wichmann-Kroll potential multiplied by r(i).
c
c       Common block parameters:
c       cl     ..... Speed of light
c       z      ..... Nuclear charge
c       ii     ..... Number of grid points
c       h, al, bt ..... Parameter of the radial grid
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_const/cl,z
        real*8 r(maxii),v_wk(maxii)
        integer ii
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Calculation Wichmann-Kroll potential times on r
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call SHAB_wk_init
c
        ii=maxii-20
	Do i=1,ii
          x=r(i)*cl
          v_wk(i)=SHAB_u_wk(x)*r(i)
	End do
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        imax=ii
        do i=ii,1,-1
	  If (dabs(v_wk(i)/r(i)).gt.1.D-40) go to 200
          imax=i
          v_wk(i)=0.d0
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
 200    v_wk(ii+3)=imax
        v_wk(ii+4)=0.d0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c        open (unit=15,file='wk.dat',status='unknown')
c        write(15,'(f22.14,16x,a)') z,'z'
c        write(15,'(i22,16x,a)') ii,'ii'
c        write(15,'(d22.14,16x,a)')  h,'h'
c        write(15,'(d22.14,16x,a)') al,'al'
c        write(15,'(d22.14,16x,a)') bt,'bt'
c        write(15,*)
c        write(15,'(17x,a,23x,a,19x)') 'r','wk'
c        do i=1,ii
c          write (15,'(i6,3d24.14)') i,r(i),v_wk(i)
c        enddo
c        write(15,*)
c        write(15,'(d22.14,16x,a)') v_wk(ii+3),'wk(ii+3)'
c        write(15,'(d22.14,16x,a)') v_wk(ii+4),'wk(ii+4)'
c        close(unit=15)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        subroutine SHAB_wk_init
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       Initialization
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        real*8 f,p,q,a
        integer s
        common /SHAB_wk_coeff/f(15,5),p(0:35),q(0:15),
     &  a(4,0:7,4),s(4,4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do m=1,15
        do n=1,5
          f(m,n) = 0.d0
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        p(0)=  2.094 001 968 517 539 d-02
        p(1)= -6.328 682 056 509 706 d-02
        p(2)=  1.656 675 047 469 226 d-01
        p(3)=  6.218 254 052 794 603 d-02
        p(4)=  8.561 450 880 149 488 d-01
        p(5)= -1.324 099 530 525 646 d+00
        p(6)=  1.022 894 990 055 742 d-01
        p(7)=  1.818 708 880 809 046 d-01
        p(8)=  5.002 561 255 243 687 d-04
        p(9)=  2.167 778 174 657 628 d-01
        p(10)= 5.453 835 857 080 859 d-01
        p(11)= 3.137 663 571 113 079 d-01
        p(12)=-1.619 774 880 547 481 d-01
        p(13)=-6.749 903 602 516 955 d-02
        p(14)=-2.935 144 138 820 853 d-04
        p(15)= 9.755 800 734 714 044 d-02
        p(16)=-1.994 016 822 313 688 d-01
        p(17)=-4.777 623 414 220 549 d-02
        p(18)= 7.503 472 520 308 676 d-03
        p(19)= 3.440 742 164 594 281 d-05

        p(20)= 2.829 420 9807 d-03
        p(21)= 1.419 551 8237 d-04
        p(22)= 1.348 983 3978 d-05
        p(23)= 2.013 400 8556 d-07
        p(24)= 8.602 918 4536 d-06
        p(25)=-1.787 016 5724 d-05
        p(26)= 2.326 960 7856 d-05
        p(27)=-1.779 846 3760 d-05
        p(28)= 8.727 474 4996 d-06
        p(29)=-2.923 560 6110 d-06
        p(30)= 6.892 936 0802 d-07
        p(31)=-1.149 169 9126 d-07
        p(32)= 1.330 283 5807 d-08
        p(33)=-1.019 098 8753 d-09
        p(34)= 4.650 804 1566 d-11
        p(35)=-9.578 592 4366 d-13
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        q(0)  = 2.09400223235d-02
        q(1)  = 1.42425723499d-02
        q(2)  = 0.76419630941d-02
        q(3)  = 0.36124746165d-02
        q(4)  = 0.15118743043d-02
        q(5)  = 0.05247168246d-02
        q(6)  = 0.01797758276d-02
        q(7)  =-0.04761316284d-02
        q(8)  = 0.19713924639d-02
        q(9)  =-0.62421312791d-02
        q(10) = 1.29259217859d-02
        q(11) =-1.89507086023d-02
        q(12) = 1.89089237003d-02
        q(13) =-1.24611625668d-02
        q(14) = 0.48793792307d-02
        q(15) =-0.08858438257d-02
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        f(5,1) = 0.3812769d+00
        f(5,2) = 0.0266795d+00
        f(5,3) = 0.0000860d+00

        f(6,1) = 0.4946161d+00
        f(6,2) = 0.0458093d+00
        f(6,3) = 0.0002918d+00
        f(6,4) = 6.d-08

        f(7,1) = 0.6150374d+00
        f(7,2) = 0.0707397d+00
        f(7,3) = 0.0007216d+00
        f(7,4) = 5.d-07

        f(8,1) = 0.7412290d+00
        f(8,2) = 0.1014082d+00
        f(8,3) = 0.0014705d+00
        f(8,4) = 23.0d-07
        f(8,5) = 1.0d-10

        f(9,1) =0.8721898d+00
        f(9,2) =0.1376368d+00
        f(9,3) =0.0026288d+00
        f(9,4) =69.d-07
        f(9,5) =2.0d-09

        f(10,1) =1.0071432d+00
        f(10,2) =0.1791838d+00
        f(10,3) =0.0042776d+00
        f(10,4) =167.0d-07
        f(10,5) =8.0d-09

        f(11,1) = 1.1454776d+00
        f(11,2) = 0.2557770d+00
        f(11,3) = 0.0064866d+00
        f(11,4) = 346.0d-07
        f(11,5) = 3.0d-08


        f(12,1) = 1.2867042d+00
        f(12,2) = 0.2771342d+00
        f(12,3) = 0.0093134d+00
        f(12,4) = 643.0d-07
        f(12,5) = 9.0d-08

        f(13,1) = 1.4304276d+00
        f(13,1) = 0.3329753d+00
        f(13,1) = 0.0128043d+00
        f(13,1) = 1097.0d-07
        f(13,1) = 2.0d-07

        f(14,1) = 1.5763244d+00
        f(14,1) = 0.3930296d+00
        f(14,1) = 0.0169953d+00
        f(14,1) = 1750.0d-07
        f(14,1) = 5.0d-07

        f(15,1) = 1.7241274d+00
        f(15,1) = 0.4570401d+00
        f(15,1) = 0.0219132d+00
        f(15,1) = 2644.0d-07
        f(15,1) = 9.0d-07
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        do i=1,4
        do m=0,7
        do n=1,4
          a(i,m,n)=0.d0
        enddo
        enddo
        enddo
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a(1,0,1) = 2.9675d+00
        a(1,0,2) = 7.8894d-01
        a(1,0,3) =-3.0201d-04
        a(1,0,4) = 4.9725d-04

        a(1,1,1) =-4.4288d+00
        a(1,1,2) =-2.9946d+01
        a(1,1,3) =-9.0441d-01
        a(1,1,4) = 2.4826d+00

        a(1,2,1) = 5.9339d01
        a(1,2,2) = 8.4881d+02
        a(1,2,3) = 5.8388d+01
        a(1,2,4) =-1.0745d+02

        a(1,3,1) =-5.6820d+02
        a(1,3,2) =-1.4888d+04
        a(1,3,3) =-1.5899d+03
        a(1,3,4) = 2.7128d+03

        a(1,4,1) = 4.6131d+03
        a(1,4,2) = 1.6003d+05
        a(1,4,3) = 2.2880d+04
        a(1,4,4) =-3.8249d+04

        a(1,5,1) =-2.6282d+04
        a(1,5,2) =-1.0259d+06
        a(1,5,3) =-1.8036d+05
        a(1,5,4) = 3.0100d+05

        a(1,6,1) = 8.9448d+04
        a(1,6,2) = 3.5980d+06
        a(1,6,3) = 7.3716d+05
        a(1,6,4) =-1.2368d+06

        a(1,7,1) =-1.3412d+05
        a(1,7,2) =-5.3067d+06
        a(1,7,3) =-1.2228d+06
        a(1,7,4) = 2.0678d+06
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a(2,0,1) =-6.64673d-01
        a(2,0,2) = 4.68647d-05
        a(2,0,3) = 2.56839d-02

        a(2,1,1) = 5.70540d-01
        a(2,1,2) =-9.83836d-04
        a(2,1,3) =-3.21021d-03

        a(2,2,1) =-9.43955d-01
        a(2,2,2) = 8.28746d-03
        a(2,2,3) = 1.01163d-02

        a(2,3,1) = 1.18405d+00
        a(2,3,2) = 5.94372d-02
        a(2,3,3) =-1.62093d-02

        a(2,4,1) =-8.91689d-01
        a(2,4,2) =-5.31431d-02
        a(2,4,3) = 1.17762d-02

        a(2,5,1) = 3.64707d-01
        a(2,5,2) = 7.52653d-03
        a(2,5,3) =-4.69814d-03

        a(2,6,1) =-6.17902d-02
        a(2,6,2) =-4.54081d-04
        a(2,6,3) = 7.76469d-04
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a(3,0,1) = 1.93743d-02
        a(3,0,2) = 1.64381d-01
        a(3,0,3) = 2.53250d-02

        a(3,1,1) =-6.43399d-01
        a(3,1,2) =-3.05149d-01
        a(3,1,3) = 3.56651d-03

        a(3,2,1) = 2.56565d-01
        a(3,2,2) = 2.31560d-01
        a(3,2,3) =-8.14980d-03

        a(3,3,1) =-1.09235d-01
        a(3,3,2) =-9.49948d-02
        a(3,3,3) = 6.21168d-03

        a(3,4,1) = 4.14226d-02
        a(3,4,2) = 2.22450d-02
        a(3,4,3) =-3.29527d-03

        a(3,5,1) =-7.28744d-03
        a(3,5,2) =-2.80168d-03
        a(3,5,3) = 7.60885d-04

        a(3,6,1) = 4.54914d-04
        a(3,6,2) = 1.47382d-04
        a(3,6,3) =-6.03993d-05
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        a(4,0,1) =-2.4734758d+09
        a(4,0,2) =-1.5115262d+04

        a(4,1,1) = 3.2270784d+09
        a(4,1,2) = 1.4945586d+04

        a(4,2,1) =-1.7207915d+09
        a(4,2,2) =-5.8361389d+03

        a(4,3,1) = 4.7635740d+08
        a(4,3,2) = 1.1476935d+03

        a(4,4,1) =-7.1414759d+07
        a(4,4,2) =-1.2098290d+02

        a(4,5,1) = 5.4305859d+06
        a(4,5,2) = 6.5439613d+00

        a(4,6,1) =-1.6479928d+05
        a(4,6,2) =-1.4296421d-01
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        s(1,1) = 1
        s(1,2) = 1
        s(1,3) = 0
        s(1,4) = 0

        s(2,1) = 0
        s(2,2) =-3
        s(2,3) = 0
        s(2,4) = 0

        s(3,1) =-1
        s(3,2) = 2
        s(3,3) = 0
        s(3,4) = 0

        s(4,1) =-13
        s(4,2) =-6
        s(4,3) = 0
        s(4,4) = 0
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
c       =================================================
        real*8 function SHAB_u_wk(x)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit real*8 (a-h,o-z)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
c       A G Fainshtein, N L Manakov and A A Nekipelov,
c       J.Phys.E: At.Mol.Opt.Phys. v.23 (1990) 559-569.
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        common /SHAB_const/cl,z
        real*8 break (5)
        parameter (pi=3.1415926535 8979323846d0)
        real*8 f,p,q,a
        integer s
        common /SHAB_wk_coeff/f(15,5),p(0:35),q(0:15),
     &  a(4,0:7,4),s(4,4)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        alz=z/cl
        dlam=1.d0-dsqrt(1.d0-alz)*dsqrt(1.d0+alz)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        if (x.lt.10.d0) then
          i=1
          if (x.gt.4.d0) then
            i=4
          else
            if (x.gt.1.5d0) then
              i=3
            else
              if (x.gt.0.15d0) i=2
            endif
          endif
c
          v3=0.d0
          if (x.le.4.25d0) then
            do n=0,8
              v3=v3+p(n)*x**(n-1)
            enddo
            do n=0,5
              v3=v3+p(n+9)*x**(n+2)*dlog(x)
            enddo
            do n=0,4
              v3=v3+p(n+15)*x**(n+3)*dlog(x)**2
            enddo
          else
            do n=0,15
              v3=v3+p(n+20)*(1.d1/x)**(2*n)
            enddo
            v3=v3/x**5
          endif
          v3=v3/cl*alz**3
c
          SHAB_u_wk=v3
c
          if (x.lt.0.15d0) then
            q5=0.d0
            do n=1,15
              q5=q5+q(n)*dlam**n
            enddo
            q5=q5*alz**3/cl
            SHAB_u_wk=SHAB_u_wk+q5/x
          endif
c
          sum=1.d0
          do n=1,4
          do m=0,7
            sum=sum+a(i,m,n)*(dlam**n)*x**(m+s(i,n))
          enddo
          enddo
          SHAB_u_wk=SHAB_u_wk/sum
c
        else
c
          SHAB_u_wk=2.d0/225.d0/x**5+59.d0/1323.d0/x**7+(1977.d0+
     &    20.d0*alz*alz)/4725.d0/x**9+(1258680.d0+34144.d0*alz*alz)/
     &    190575/x**11+(1960420032.d0+93618070.d0*alz*alz+
     &    96096.d0*alz**4)/12297285.d0/x**13
          SHAB_u_wk=SHAB_u_wk*alz**3
c
          prod=576.d0/x**13
          do m=5,15
            prod=prod*m*m/(x*x)
            sum=0.d0
            do n=1,5
              sum=sum+f(m,n)*alz**(2*n+1)
            enddo
            SHAB_u_wk=SHAB_u_wk+sum*prod
          enddo
          SHAB_u_wk=SHAB_u_wk/(pi*cl)
        endif
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
