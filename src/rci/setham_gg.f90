!>
!!
!! Comment from the original `SETHAM_GG` header:
!!
!! > Sets up the Hamiltonian matrix and determines the average energy.  *
!! >
!! > Serial I/O moved out; able to run on single/multiple processors
!! > For this purpose a new common /setham_to_genmat2/ is created                                                           *
      SUBROUTINE SETHAM (myid, nprocs, jblock, ELSTO,ICSTRT, nelmntt)
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: NNNP, KEYORB
      USE memory_man
!cjb iso_fortran_env :: output_unit, error_unit
      use, intrinsic :: iso_fortran_env
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE bcore_C, ONLY: icore
      USE bilst_C
      USE buffer_C
      USE coeils_C
      USE decide_C
      USE debug_C, only: IBUG1
      USE def_C
      USE eigv_C
      USE iccu_C
      USE grid_C
      USE hmat_C
      USE keilst_C
      USE ncdist_C
      USE orb_C, ONLY: ncf,  nw
      USE prnt_C
      USE stat_C
      USE stor_C
      USE tatb_C
      USE vinlst_C
      USE vpilst_C
      USE wave_C
      USE cteilsrk_C
      USE blim_c
      USE where_C
      USE eigvec1_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      USE ichop_I
      use cord_I
      use grasp_rciqed_qed, only: slfint
      use grasp_cimatrixelements
      ! Matrix elements smaller than cimatrixelements%cutoff are not accumulated
      IMPLICIT NONE
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER(LONG)       :: NELMNTT
      INTEGER             :: JBLOCK, ICSTRT
      INTEGER, INTENT(IN) :: MYID, NPROCS
      REAL(DOUBLE)        :: ELSTO
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      REAL(DOUBLE), DIMENSION(NNNW) :: tshell
      REAL(DOUBLE) :: breit_core, breit_noncore

      INTEGER, PARAMETER :: KEY = KEYORB

      INTEGER :: ipi, ipj, inc1, inc2, kt, ipt, incor, ncoec, nctec, &
                 i, j, nmcbp, ncore, ic, nelc, irstart, ir, ia, ib,  &
                 itype, nctei, iia
       REAL(DOUBLE) :: elemnt, precoeff, tcoeff, vcoeff, contr, t1, t2, t3
!-----------------------------------------------------------------------
      nelmnt = nelmntt

      IF (IPRERUN .EQ. 2) THEN
         DO IPI = 1,NVEC
            DO IPJ = 1,NCF
               WRITE (*,*) IPI,IPJ,EVEC1(IPJ+(IPI-1)*NCF)
            ENDDO
         ENDDO
      ENDIF
!
!   Allocate storage to arrays in COMMON/BUFFER/; these are
!   used for the Coulomb and transverse two-electron integrals
!
      CALL ALCBUF (1)

!     ...Locals
      CALL alloc (emt, ncf, 'EMT','SETHAM' )
      CALL alloc (irow, ncf, 'IROW', 'SETHAM')
!
      INC1 = 1
      INC2 = 1
!
!   Initialisations for contributions from the Dirac-Coulomb
!   operator
!
      KT  = 0
      IPT = 1
!
      INCOR = 1

      NCOEC = 0
!
      NCTEC   = 0

      IF (LTRANS) THEN

!        ...Initialisations for transverse interaction correction
         DO 2 I = 1, NW
            ICORE(I) = 0
            DO J = 1, NCF
               IF (ICHOP (I,J) .LE. 0) GOTO 2
            ENDDO
            ICORE(I) = 1
    2    CONTINUE

         NMCBP = 0
         NCORE = 0
      ENDIF

! Loop over rows of the Hamiltonian matrix - distributed

      DO 10 ic = icstrt, ncf, nprocs

         NELC = 0    ! counter - Number of non-zeros of this row

! Loop over columns of the current row

         irstart = 1
         DO 85 IR = irstart, IC

! PER
            IF (LFORDR .AND. (IR .GT. ICCUT(1))) THEN
               IF (IR.NE.IC) CYCLE
            END IF
! PER

            ELEMNT = 0.D0     ! accumulates various contributions to H

!
!   Generate the integral list for the matrix element of the
!   one-body operators
!
            IF (IPRERUN .EQ. 1) THEN
               INC1 = 0
               INC2 = 0
               IF (IC.LE.NCSFPRE .OR. IC.EQ.IR) THEN
                  INC1 = 1
               ENDIF
            ENDIF

            IF (IPRERUN .EQ. 2) THEN
!
!   Diagonal elements are always included
!   Off diagonal elements are included only if the
!   products of the weights from the prerun are larger
!   than the cutoff.
!
               IF (IC .EQ. IR) THEN
                  INC1 = 1
                  INC2 = 1
               ELSE
                  INC1 = 0
                  INC2 = 0
               ENDIF
               DO IPI = 1,NVEC
                  PRECOEFF =                                           &
                   DABS(EVEC1(IC+(IPI-1)*NCF)*EVEC1(IR+(IPI-1)*NCF))
                  IF (PRECOEFF .GT. COEFFCUT1) INC1 = 1
                  IF (PRECOEFF .GT. COEFFCUT2) INC2 = 1
               ENDDO
            ENDIF

!            ...INC1.EQ.1 ------------>
         IF (INC1 .EQ. 1) THEN   !inc1 is always 1 without PRE-RUN

             !   Accumulate the contribution from the one-body operators:
             !   kinetic energy, electron-nucleus interaction; update the
             !   angular integral counter
           CALL ONESCALAR(IC,IR,IA,IB,TSHELL)
           ELEMNT = ELEMNT + dirac_potential(IC, IR, IA, IB, TSHELL)
           IF(LNMS) ELEMNT = ELEMNT + nms(IC, IR, IA, IB, TSHELL)
           IF(LVP) ELEMNT = ELEMNT + qed_vp(IC, IR, IA, IB, TSHELL)

           IF (IA .NE. 0) THEN
              IF (IA == IB) THEN
                  DO IA = 1, NW
                      IF(ABS (TSHELL(IA)) .GT. CUTOFF) NCOEC = NCOEC + 1
                  ENDDO
              ELSE IF(ABS(TSHELL(1)) .GT. CUTOFF) THEN
                 NCOEC = NCOEC + 1
              ENDIF
           ENDIF
!
!
!
         IBUG1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
         NVCOEF = 0
         CALL RKCO_GG (IC, IR, CORD, INCOR, 1)
         ELEMNT = ELEMNT + coulomb_cached(ic, ir)
         IF(LSMS) ELEMNT = ELEMNT + sms_cached(ic, ir)

         DO I = 1, NVCOEF
            VCOEFF = DBLE(COEFF(I))
            IF (DABS (VCOEFF) .GT. CUTOFF) THEN
               NCTEC = NCTEC + 1
           ENDIF
         END DO
!
         IBUG1 = 0

         ENDIF  !inc1 is always 1 without PRE-RUN
!            ...INC1.EQ.1 <------------
!***********************************************************************
!            ...LTRANS .AND. (INC2.EQ.1) ------------>
         IF (LTRANS .AND. (INC2.EQ.1)) THEN
            !IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
            call breit_split(ic, ir, breit_core, breit_noncore)
            ELEMNT = ELEMNT + breit_noncore
            ELSTO = ELSTO + breit_core
            ! ELSTO is a constant over all diagonals, thus its contribution to
            ! the total energy can be added later
            ! IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO

            DO I = 1, NVCOEF
               IF (DABS (COEFF(I)) > CUTOFF) THEN
                  NMCBP = NMCBP + 1
                  IF(.not.(LABEL(6,I) > 0)) then
                      ! It comes here only when ic=ir=1
                      ! clue: rkco<-breid<-talk<-label(6,i)
                      NCORE = NCORE + 1
                  endif
               ENDIF
            ENDDO
            IBUG1 = 0
         ENDIF
!            ...LTRANS .AND. (INC2.EQ.1) <------------
!***********************************************************************
!
! Store this element if it is diagonal or its value is greater than
! CUTOFF
!
         IF ( (IR == IC) .OR. (DABS (ELEMNT) > CUTOFF) ) THEN
            NELC       = NELC + 1
            EMT(NELC)  = ELEMNT
            IROW(NELC) = IR
         ENDIF
!
   85    CONTINUE

         ! If LSE, then add self-energy to the diagonal of the CI matrix.
         IF(LSE) EMT(NELC) = EMT(NELC) + qed_se_mohr(IC, slfint)

!   This column is done; write it to disk
!
         WRITE (imcdf) NELC, ELSTO, (EMT(IR), IR = 1, NELC),       &
                                   (IROW(IR), IR = 1, NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This EAV (and the above EMT) does not have ELSTO.
         EAV = EAV + EMT(NELC)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Output IC on screen to show how far the calculation has proceeded
!
!-----------------------------------------------------------------------
!cjb      if (mod(ic-1,100).eq.0) then
!cjb         PRINT '(A4,I10,A9,I10,A9,I5,A8,I5)', 'Row ', ic,       &
!cjb                 ' nnonz = ', nelc,                             &
!cjb                 ' block = ', jblock,' myid = ', myid
!cjb      ENDIF
!
!cjb          CALL CONVRT (IC,CNUM,LCNUM)
!cjb          if (mod(IC,100).eq.0) then
!cjb             PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
!cjb          end if
!cjb
!
!cjb MPI progress begin --------------------------------------------------
!
!cjb just started
       IF (IC .LE. nprocs) then
            PRINT '(A6,I10,A9,I10,A9,I5,A8,I5)', 'Start ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
!
!cjb progress
       IF (IC .GT. nprocs .and. IC .LT. NCF-nprocs .and. &
           MOD (IC-1,100*nprocs) .EQ. myid) then
         if (myid .eq. 0) then
            PRINT '(A5,I10,A9,I10,A9,I5,A8,I5)', 'Done ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         else
            PRINT '(A5,I10,A9,I10,A9,I5,A8,I5)', 'Row  ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
         endif
       endif
!
!cjb almost finished
       IF (IC .GT. NCF-nprocs) then

            PRINT '(A7,I10,A9,I10,A9,I5,A8,I5)', 'Finish ', ic,       &
                    ' nnonz = ', nelc,                             &
                    ' block = ', jblock,' myid = ', myid
        flush output_unit
        flush error_unit
       endif
!cjb MPI progress end ----------------------------------------------------
!
!
!   Update the counter for the total number of elements
!
         NELMNT = NELMNT + NELC
!
   10 CONTINUE
!***********************************************************************
!
!   Deallocate storage for the arrays in /BUFFER/
!
      CALL ALCBUF (3)

!     ...Locals
      CALL DALLOC (EMT, 'EMT', 'SETHAM')
      CALL DALLOC (IROW, 'IROW', 'SETHAM')

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      CUTOFFtmp = CUTOFF
      NCOEItmp = NCOEI
      NCOECtmp = NCOEC
      NCTEItmp = NCTEI
      NCTECtmp = NCTEC
      NTPItmp = NTPI
      NMCBPtmp = NMCBP
      NCOREtmp = NCORE
      NVPItmp = NVPI
      NKEItmp = NKEI
      NVINTItmp = NVINTI
      NELMNTtmp = NELMNT
      NCFtmp = NCF

      RETURN
      END SUBROUTINE SETHAM
