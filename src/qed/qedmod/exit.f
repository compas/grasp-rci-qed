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
! Originally from: f/potgen/potgen_main.f, f/potuse/potuse1.f, f/potgen/potuse2.f
!
c       =================================================
        subroutine SHAB_exit1
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        call exit(1)
c       - - - - - - - - - - - - - - - - - - - - - - - - -
        return
        end
