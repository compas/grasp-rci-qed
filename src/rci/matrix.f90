!***********************************************************************
!                                                                      *
      SUBROUTINE MATRIX (ncore, j2max)
!                                                                      *
!   This SUBROUTINE calls routines to  form the  Hamiltonian  matrix   *
!   and to diagonalise it.   The total angular momenta and parity of   *
!   the ASF are found;  the eigenvectors are chosen so that the sign   *
!   of the largest element is positive.                                *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, HMOUT, IQ              *
!                        WGHTD5.                                       *
!               [RCI92]: MANEIG, QED, SETHAM.                          *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Modified by Xinghong He               Last revision: 31 Dec 1997   *
!   Block version Xinghong He             Last revision: 12 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE parameter_def,   ONLY: NNNW
      USE memory_man
      USE decide_C
      USE DEBUG_C
      USE def_C
      USE eigv_C
      USE hblock_C
      USE hmat_C
      USE iounit_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE orb_C, ONLY: ncf, nw, iqa
      USE prnt_C
      USE stat_C
      USE wave_C
      USE blim_C
      USE eigvec1_C
      USE iccu_C
      USE cteilsrk_C
      USE coeils_C
      USE bilst_C
      USE keilst_C
      USE vinlst_C
      USE fposition_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE auxblk_I
      USE lodcslmpi_I
      USE qed_slfen_I
      USE genmat_I
      USE genmat2_I
      USE hmout_I
      USE iq_I
      USE maneig_I
      USE engout_I
      USE wghtd5_I
      USE mpi_C

      use grasp_rciqed, only: IMCDF => res_unit, setype, sematrix
      use grasp_rciqed_qed, only: qedse
      use grasp_rciqed_cimatrixelements, only: qed_se_mohr

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: ncore
      INTEGER,  INTENT(IN):: j2max
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=8) :: CNUM
      REAL(DOUBLE) :: elsto, eau, ecm, eev, elemnt
      REAL(DOUBLE), DIMENSION(:), pointer :: etot
      INTEGER(LONG) :: nelmnt_a
      INTEGER :: iiatjpo, iiaspar
      INTEGER :: i, j, irestart, ncminpas, jblock, ic, ip, mode, k
!-----------------------------------------------
!     ...Common to all blocks - place here to save CPU time
      CALL auxblk (j2max)
      ! Populate the global array in grasp_rciqed with the matrix elements of
      ! the self-energy operator picked in GETCID.
      if(LSE) then
          call qedse(setype, sematrix)
      endif

!***************************************************************
!      Loop over blocks
!***************************************************************
      ncminpas = 0

      DO 100 jblock = 1, nblock
         ncf    =   ncfblk(jblock)
         nvec   =   nevblk(jblock)
         nvecmx = ncmaxblk(jblock)
         iccut(1)  = iccutblk(jblock)
         !.. Determine position of the previous block in the .res file

         nposition = 7 + nw + nw  ! File position of the previous block
                                  ! in the .res file
         DO i = 1, jblock - 1
            j = (ncfblk(i) - myid - 1 + nprocs) / nprocs
            IF (ncfblk(i) .LT. nprocs) j = ncfblk(i) / (myid+1)
            nposition = nposition + j + 1
         ENDDO

         !.. SETHAM does not handle this extrem case
         IF (nprocs .GT. NCF)                                          &
                         CALL stopmpi ('matrix: too many nodes', myid)

!        ...Obtain ivec() from iccmin()
         IF (nvec .GT. 0) THEN
            CALL alloc (ivec, nvec, 'IVEC', 'MATRIX')
            DO i = 1, nvec
               ivec(i) = iccmin(i+ncminpas)
            ENDDO
            ncminpas = ncminpas + nvec
         ENDIF

!        ...These 3 were allocated in lodcsh2 and deallocated at the end
!        ... of this routine and in the setham. In this block version,
!        ... both allocation and deallocation are placed here. See the
!        ... following goto 80 for reason.
         CALL ALLOC (IQA, NNNW, ncf, 'IQA', 'MATRIX')
         CALL ALLOC (JQSA, NNNW, 3, ncf, 'JQSA', 'MATRIX')
         CALL ALLOC (JCUPA, NNNW, ncf, 'JCUPA', 'MATRIX')

!      ...Load CSF list of the current block
        CALL lodcslmpi (21, ncore, jblock)

         IF (nvec <= 0) THEN
            IF (myid .EQ. 0) WRITE (25) jblock, ncf, nvec, 999, 999
!           ...Generate H matrix of the current block and write to disk
!           ...eav, nelmnt are also obtained from genmat
            CALL genmat (jblock, myid, nprocs, elsto, irestart)
             ! This call is optional, meaning the matrix of this block
               ! does not need to be computed. But don't comment it out
               ! since other programs may assume thet existence of them.
            CALL genmat2 (irestart, nelmnt_a, elsto)
            GOTO 80 ! need to clear memory
         ENDIF
!        ------------------------
         CALL genmat (jblock, myid, nprocs, elsto, irestart)
         CALL genmat2 (irestart, nelmnt_a, elsto)
!
!   Allocate and deallocate memory for the mixing coefficients
!   from the prerun
!
      IF (IPRERUN.EQ.1) CALL ALLOC (EVEC1,NCF*NVEC, 'EVEC1', 'MATRIX')
      IF (IPRERUN.EQ.2) CALL DALLOC (EVEC1, 'EVEC1', 'MATRIX')

! Since maneig needs both nelmnt and nelmnt_a, # of non-zero
! matrix elements computed by current node and the total #
! from all nodes, a new var nelmnt_a is created rathan than
! using the one in the common block.

!        write (*,*) 'Matrix61,myid=',myid
        CALL MANEIG (iiatjpo, iiaspar, nelmnt_a)
!GG        CALL MANEIG (iiatjpo, iiaspar, nelmnt_a, imethod)
!GG         CALL MANEIG (iiatjpo, iiaspar)
!
!   Write out eigenvalues (ENGOUT), dominant components of the
!   eigenvectors (WGHTD5) to stream 24; write out ASF symmetries,
!   eigenvalues  and eigenvectors to RCI92 mixing coefficients file.
!   EAV and ELSTO are added back to energy here
!
        IF (myid .EQ. 0) THEN
! ELSTO has never been in Hamiltonian matrix, yet it was
! added to EAV which was later substracted from H. Thus at
! this point, EAV is correct (it has ELSTO added), EVAL()
! need ELSTO and the correct EAV.
           IF (NCF > 1) then
              DO i = 1, NVEC
                 EVAL(i) = EVAL(i) + ELSTO
              ENDDO
           END IF

           CALL ENGOUT (EAV,EVAL,IiATJPO,iIASPAR,IVEC,NVEC,3)
           CALL WGHTD5 (iiatjpo, iiaspar)

!          ...Write ASF symmetries, eigenvalues, and eigenvectors to RCI92
!          ...MIXing coefficients File; close the file; print a report
           WRITE (25) jblock, ncf, nvec, iiatjpo, iiaspar
           WRITE (25) (ivec(i), i = 1, nvec)
           WRITE (25) EAV,(EVAL(I),I = 1,NVEC)
           WRITE (25) ((EVEC(I+(J-1)*NCF),I = 1,NCF),J = 1,NVEC)

           PRINT *, 'RCI90 MIXing coefficients File generated.'
        ENDIF
!
!   Save the mixing coefficients from the prerun
!
      IF (IPRERUN .EQ. 1) THEN
         DO J = 1, NVEC
            DO I = 1, NCF
               EVEC1(I+(J-1)*NCF) = EVEC(I+(J-1)*NCF)
            ENDDO
         ENDDO
      ENDIF

!        ...Locals
         CALL dalloc (ivec, 'IVEC', 'MATRIX')
!        ...Allocated in maneig
         CALL dalloc (eval, 'EVAL', 'MATRIX')
         CALL dalloc (evec, 'EVEC', 'MATRIX')

   80    CONTINUE

!        ...Locals
         CALL dalloc (IQA, 'IQA', 'MATRIX')
         CALL dalloc (JQSA, 'JQSA', 'MATRIX')
         CALL dalloc (JCUPA, 'JCUPA', 'MATRIX')

  100 CONTINUE
!
!   Close the restart files; nothing will be added to them now
!
      CLOSE (imcdf)
!     ...Clean up
      CALL dalloc (ncfblk, 'NCFBLK', 'MATRIX')
      CALL dalloc (nevblk, 'NEVBLK', 'MATRIX')
      CALL dalloc (ncmaxblk, 'NCMAXBLK', 'MATRIX')
      CALL dalloc (iccutblk, 'ICUTTBLK', 'MATRIX')
      CALL dalloc (iccmin, 'ICCMIN', 'MATRIX') ! allocated in items as pnccmin

      CALL dalloc (VALTEIRK, 'VALTEIRK', 'MATRIX') ! allocated in genintrk
      CALL dalloc (INDTEIRK, 'INDTEIRK', 'MATRIX') ! allocated in genintrk
!
!   Deallocate storage for the integral lists from the
!   Dirac-Coulomb operator; the storage was allocated
!   in IABINT and RKINTC
!
      IF (NCOEI .GT. 0) THEN
         CALL DALLOC (INDOEI, 'INDOEI', 'MATRIX')
         CALL DALLOC (VALOEI, 'VALOEI', 'MATRIX')
      ENDIF
!
!   Deallocate storage for the integral lists from the
!   transverse photon interaction operator; this storage
!   was allocated in BRINT1, brint2,...brint6
!
      IF (LTRANS) THEN
         IF (NTPI(1) .GT. 0) THEN
            CALL DALLOC (INDTP1, 'INDTP1', 'MATRIX')
            CALL DALLOC (VALTP1, 'VALTP1', 'MATRIX')
         ENDIF
         IF (NTPI(2) .GT. 0) THEN
            CALL DALLOC (INDTP2, 'INDTP2', 'MATRIX')
            CALL DALLOC (VALTP2, 'VALTP2', 'MATRIX')
         ENDIF
         IF (NTPI(3) .GT. 0) THEN
            CALL DALLOC (INDTP3, 'INDTP3', 'MATRIX')
            CALL DALLOC (VALTP3, 'VALTP3', 'MATRIX')
         ENDIF
         IF (NTPI(4) .GT. 0) THEN
            CALL DALLOC (INDTP4, 'INDTP4', 'MATRIX')
            CALL DALLOC (VALTP4, 'VALTP4', 'MATRIX')
         ENDIF
         IF (NTPI(5) .GT. 0) THEN
            CALL DALLOC (INDTP5, 'INDTP5', 'MATRIX')
            CALL DALLOC (VALTP5, 'VALTP5', 'MATRIX')
         ENDIF
         IF (NTPI(6) .GT. 0) THEN
            CALL DALLOC (INDTP6, 'INDTP6', 'MATRIX')
            CALL DALLOC (VALTP6, 'VALTP6', 'MATRIX')
         ENDIF
      ENDIF
!
!   Deallocate storage for the nuclear motional energy integral
!   lists; this was allocated in KEINT and VINT
!
      IF (LNMS) THEN
         IF (NKEI .GT. 0) THEN
            CALL DALLOC (INDKEI, 'INDKEI', 'MATRIX')
            CALL DALLOC (VALKEI, 'VALKEI', 'MATRIX')
         ENDIF
      ENDIF
      IF (LSMS) THEN
         IF (NVINTI .GT. 0) THEN
            CALL DALLOC (INDTEI, 'INDTEI', 'MATRIX')
            CALL DALLOC (VALTEI, 'VALTIE', 'MATRIX')
         ENDIF
      ENDIF

      CALL dalloc (PF, 'PF', 'MATRIX') ! lodrwf or lodres
      CALL dalloc (QF, 'QF', 'MATRIX') ! lodrwf or lodres

      RETURN

      END
