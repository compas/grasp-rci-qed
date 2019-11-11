!***********************************************************************
!                                                                      *
      REAL(KIND(0.0D0)) FUNCTION RATDEN (P, Q, MTPO, NP, KAPPA, Z)
!                                                                      *
!   This function calculates the ratio of occupation of the orbital in *
!   the P and Q arrays and the hydrogenic orbital defined by NP, KAPPA *
!   and Z in the radial range 0 < r < 0.0219 a.u.                      *
!                                                                      *
!   I.e. in other words, it integrates the densities of the two        *
!   orbitals in that range and the returns their ratio.                *
!                                                                      *
!   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Switches:
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: NNNP
      USE grid_C
      USE tatb_C, ONLY: mtp, ta
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dcbsrw_I
      USE quad_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MTPO
      INTEGER  :: NP
      INTEGER  :: KAPPA
      REAL(DOUBLE)  :: Z
      REAL(DOUBLE), INTENT(IN) :: P(NNNP)
      REAL(DOUBLE), INTENT(IN) :: Q(NNNP)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MTPH, I, K
      REAL(DOUBLE) :: EH, PZH, RESULT, RESULT1
      REAL(DOUBLE), DIMENSION(NNNP) :: PH, QH
!-----------------------------------------------
!
!
!
!   Set up the hydrogenic orbital
!
      CALL DCBSRW (NP, KAPPA, Z, EH, PZH, PH, QH, MTPH)
!
!   Compute the overlap
!
      MTP = MIN(MTPH,MTPO)
      DO I = 2, MTP
         IF (RP(I) > 0.0219) CYCLE
         K = I
      END DO
      MTP = K
      TA(1) = 0.0D00
      TA(2:MTP) = (P(2:MTP)*P(2:MTP)+Q(2:MTP)*Q(2:MTP))*RP(2:MTP)
!         TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
      CALL QUAD (RESULT)
      TA(2:MTP) = (PH(2:MTP)*PH(2:MTP)+QH(2:MTP)*QH(2:MTP))*RP(2:MTP)
      CALL QUAD (RESULT1)
!
      RATDEN = RESULT/RESULT1
!
      RETURN
      END FUNCTION RATDEN
