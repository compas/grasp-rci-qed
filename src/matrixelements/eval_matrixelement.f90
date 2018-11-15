! module eval_matrixelement_common
!     use g2k_parameters
!     implicit real*8 (a-h, o-z)
!
!     private
!
!     EXTERNAL CORD
!     EXTERNAL BREID
!
!     COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
!     POINTER (PLABEL,LABEL(6,1))
!     POINTER (PCOEFF,COEFF(1))
!     public :: NVCOEF, COEFF, LABEL
!
!     POINTER (PNEVAL,EVAL(1))
!     POINTER (PNEVEC,EVEC(1))
!     POINTER (PNIVEC,IVEC(1))
!     POINTER (PIATJP,IATJPO(1))
!     POINTER (PIASPA,IASPAR(1))
!
!     COMMON/DEF1/EMN,IONCTY,NELEC,Z &
!           /EIGVAL/EAV,PNEVAL &
!           /EIGVEC/PNEVEC &
!           /ORB2/NCF,NW,PNTRIQ,NOB,NCFBL(100),NEVBL(100) &
!           /PRNT/NVEC,PNIVEC,NVECMX &
!           /SYMA/PIATJP,PIASPA
!
!     public :: CORD, BREID
!     public :: NCF, NW
!     public :: NCFBL
!     public :: NVEC
!     public :: EVEC
!     public :: EMN
!
! end module eval_matrixelement_common

!> Calculates the matrix element V(i, j)
function eval_matrixelement(ic, ir, hamcache)
    use, intrinsic :: ieee_arithmetic

    use grasp_kinds, only: real64, dp
    !use g2k_parameters
    use g2k_librci, only: hamiltonian_cache, matrixelement
    !use eval_matrixelement_common
    use buffer_C
    use def_C, only: EMN
    use orb_C
    use prnt_C
    use syma_C
    use eigv_C
    use eigvec1_C

    use cord_I
    use breid_I
    use rkint_I
    use rkintc_I
    use rkco_gg_I
    use brint1_I
    use brint2_I
    use brint3_I
    use brint4_I
    use brint5_I
    use brint6_I
    use onescalar_I
    use iabint_I
    use keint_I
    use vint_I
    use vpint_I
    implicit none

    integer, intent(in) :: ic, ir
    type(hamiltonian_cache), intent(in) :: hamcache
    type(matrixelement) :: eval_matrixelement

    ! Matrix elements smaller than CUTOFF are not accumulated. This is from
    ! matrix.f in rci_mpi.
    real(real64), parameter :: cutoff = 1d-20

    integer :: incor, itype
    real(real64) :: tshell(NW), ATWINV

    integer :: i, j, k, l, m, setype
    integer m_local
    real(real64) :: result, result1, result2

    eval_matrixelement%se(1) = ieee_value(eval_matrixelement%se(1), ieee_quiet_nan)
    eval_matrixelement%se(2) = ieee_value(eval_matrixelement%se(2), ieee_quiet_nan)
    eval_matrixelement%se(3) = ieee_value(eval_matrixelement%se(3), ieee_quiet_nan)
    eval_matrixelement%se(4) = ieee_value(eval_matrixelement%se(4), ieee_quiet_nan)

    incor = 1 ! set to this in setham_gg.f
    ATWINV = 1.0_dp/EMN ! this is calculated in setham_gg.f
    print *, EMN, ATWINV

    ! --------------------------------------------------------------------------
    ! This part calculates the one-particle matrix elements. Based on setham_gg.f
    ! from rci_mpi.
    call ONESCALAR(ic, ir, k, l, tshell)

    ! It appears that when k==0, then there is no contribution? So, following
    ! setham_gg.f, we only do something when k /= 0.\
    eval_matrixelement%diracpot = 0.0_dp
    if(k /= 0 .and. k == l) then
        ! If k==l, then it seems that we have to add just the diagonal elements
        ! of the single particle matrix elements, weighted with TSHELL, to the
        ! matrix element.
        print *, "k / l", k, l
        do m = 1, NW
            if(abs(tshell(m)) <= cutoff) cycle
            m_local = m
            call IABINT(m_local, m_local, result)
            eval_matrixelement%diracpot = eval_matrixelement%diracpot + result * tshell(m)
            print *, m, tshell(m), result
        enddo
    elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
        print *, "k / l", k, l
        call IABINT(k, l, result)
        eval_matrixelement%diracpot = eval_matrixelement%diracpot + result * tshell(1)
        print *, "tshell(1)", 1, tshell(1), result
    else
        print *, "k", k
    endif
    print '(": eval_matrixelement%diracpot -> ",d22.15)', eval_matrixelement%diracpot

    ! --------------------------------------------------------------------------
    ! QED parts of the matrix element. The loops have the same logic as for the
    ! 1-particle DC.
    !
    ! TODO: The two VP contributions should be split up.
    eval_matrixelement%vp = 0.0_dp
    if(k /= 0 .and. k == l) then
        print *, "k / l", k, l
        do m = 1, NW
            if(abs(tshell(m)) <= cutoff) cycle
            m_local = m
            call VPINT(m_local, m_local, result)
            eval_matrixelement%vp = eval_matrixelement%vp + result * tshell(m)
            print *, m, tshell(m), result
        enddo
    elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
        print *, "k / l", k, l
        call VPINT(k, l, result)
        eval_matrixelement%vp = eval_matrixelement%vp + result * tshell(1)
        print *, "tshell(1)", 1, tshell(1), result
    else
        print *, "k", k
    endif
    print '(": eval_matrixelement%vp -> ",d22.15)', eval_matrixelement%vp

    ! QED self-energy
    ! do setype = 1, 4
    !     eval_matrixelement%se(setype) = 0.0_dp
    !     if(k /= 0 .and. k == l) then
    !         print *, "k / l", k, l
    !         do m = 1, NW
    !             if(abs(tshell(m)) <= cutoff) cycle
    !             eval_matrixelement%se(setype) = eval_matrixelement%se(setype) + hamcache%sematrices(setype, m, m) * tshell(m)
    !             print *, m, tshell(m), hamcache%sematrices(setype, m, m)
    !         enddo
    !     elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
    !         print *, "k / l", k, l
    !         eval_matrixelement%se(setype) = eval_matrixelement%se(setype) + hamcache%sematrices(setype, k, l) * tshell(1)
    !         print *, "tshell(1)", 1, tshell(1), hamcache%sematrices(setype, k, l)
    !     else
    !         print *, "k", k
    !     endif
    !     print '(": eval_matrixelement%se[",i0,"] -> ",d22.15)', setype, eval_matrixelement%se(setype)
    ! enddo

    ! --------------------------------------------------------------------------
    ! Normal mass shift (NMS)
    eval_matrixelement%nms = 0.0_dp
    if(k /= 0 .and. k == l) then
        print *, "k / l", k, l
        do m = 1, NW
            if(abs(tshell(m)) <= cutoff) cycle
            m_local = m
            call KEINT(m_local, m_local, result)
            eval_matrixelement%nms = eval_matrixelement%nms + result * ATWINV * tshell(m)
            print *, m, tshell(m), result
        enddo
    elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
        print *, "k / l", k, l
        call KEINT(k, l, result)
        eval_matrixelement%nms = eval_matrixelement%nms + result * ATWINV * tshell(1)
        print *, "tshell(1)", 1, tshell(1), result
    else
        print *, "k", k
    endif
    print '(": eval_matrixelement%nms -> ",d22.15)', eval_matrixelement%nms

    ! --------------------------------------------------------------------------
    ! And here is the two particle part of DC. Still from setham_gg.f
    !
    ! The is this comment in setham_gg.f about this part of the code:
    !
    ! > Accumulate the contributions from the two-electron
    ! > Coulomb operator and the mass polarisation; the latter
    ! > is computed first because the orbital indices may be
    ! > permuted by RKINTC
    !
    ! The CORD variable appears to be a procedure from lib92.
    !
    ! Somewhere deep in the dependency tree probably, RKCO_GG needs the TALK,
    ! CXK and SNRC routines that are _defined locally_ in rangular and rci for
    ! what ever reason. They now live in rcicommon.
    NVCOEF = 0 ! this is done in SETHAM and then the variable used later
    CALL RKCO_GG(ic, ir, CORD, incor, 1)
    print '("incor = ",i0)', incor
    print '("NVCOEF = ",i0)', NVCOEF
    do i = 1, NVCOEF
        print '("> ",6i5)', i, LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i)
    enddo

    ! Coulomb interaction
    eval_matrixelement%coulomb = 0.0_dp
    do i = 1, NVCOEF
        print '("[",i0,"] eval_matrixelement%coulomb")', i
        if(abs(COEFF(i)) < cutoff) cycle
        print *, " > LABEL", LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i)
        !call RKINTC(LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i), result)

        call rkintc_safe(LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i), result)
        print *, "RKINTC", i, COEFF(i), result, result * COEFF(i)
        print *, " > LABEL", LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i)
        eval_matrixelement%coulomb = eval_matrixelement%coulomb + result * COEFF(i)
        print '(" eval_matrixelement%coulomb -> ",d22.15)', eval_matrixelement%coulomb
    enddo
    print '(": eval_matrixelement%coulomb -> ",d22.15)', eval_matrixelement%coulomb

    ! Special mass shift (SMS)
    eval_matrixelement%sms = 0.0_dp
    print '("ATWINV = ",d20.10)', ATWINV
    do i = 1, NVCOEF
        if(abs(COEFF(i)) < cutoff) cycle
        if(LABEL(5,I) == 1) then
            print '("[",i0,"] eval_matrixelement%sms")', i
            print *, " > LABEL", i, LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i)
            call VINT(LABEL(1,i), LABEL(3,i), result1)
            print '("VINT(i=",i0,"; ",i0,", ",i0,") -> ",d22.14)', i, LABEL(1,i), LABEL(3,i), result1
            call VINT(LABEL(2,i), LABEL(4,i), result2)
            print '("VINT(i=",i0,"; ",i0,", ",i0,") -> ",d22.14)', i, LABEL(2,i), LABEL(4,i), result2
            eval_matrixelement%sms = eval_matrixelement%sms - result1 * result2 * ATWINV * COEFF(i)
            print '(" eval_matrixelement%sms -> ",5d22.14)', COEFF(i), result1, result2, ATWINV, &
                result1 * result2 * ATWINV * COEFF(i)
        endif
    enddo
    print '(": eval_matrixelement%sms(",i0,", ",i0,") -> ",d22.14)', ic, ir, eval_matrixelement%sms

    ! --------------------------------------------------------------------------
    ! The Breit interaction part of the matrix element
    eval_matrixelement%breit = 0.0_dp
    NVCOEF = 0 ! needs to be reset?
    CALL RKCO_GG(ic, ir, BREID, 1, 2)
    do i = 1, NVCOEF
        if(abs(COEFF(i)) < cutoff) cycle

        itype = abs(LABEL(6,i))
        if(itype == 1) then
            call BRINT1(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        elseif(itype == 2) then
            call BRINT2(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        elseif(itype == 3) then
            call BRINT3(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        elseif(itype == 4) then
            call BRINT4(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        elseif(itype == 5) then
            call BRINT5(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        elseif(itype == 6) then
            call BRINT6(LABEL(1,I),LABEL(2,I),LABEL(3,i),LABEL(4,i),LABEL(5,i), result)
        endif
        ! Note: itype == 0 is also possible.

        print *, "BRINT*", i, COEFF(i), result, result * COEFF(i)
        print *, " > LABEL", LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i), LABEL(6,i)
        if(itype > 0) then
            eval_matrixelement%breit = eval_matrixelement%breit + result * COEFF(i)
        else
            ! Comment from the original setham_gg routine:
            !
            ! > ...It comes here only when ic=ir=1
            ! > clue: rkco<-breid<-talk<-label(6,i)
            !
            ! And there also was the following code in this else block:
            !
            !     NCORE = NCORE + 1
            !     ELSTO = ELSTO + result * COEFF(i)
            !
            ! However, this code does not look relevant. It sets up some stuff
            ! for genmat two.. but does not seem to be able to affect matrix
            ! element calculations.
        endif
        print '("[",i0,"] eval_matrixelement%breit -> ",d22.15)', i, eval_matrixelement%breit
    enddo
    print '(": eval_matrixelement%breit -> ",d22.15)', eval_matrixelement%breit

    !eval_matrixelement = eval_matrixelement%diracpot + eval_matrixelement%coulomb + eval_matrixelement%breit
    !print '("matrixelement -> ", d22.15)', matrixelement

contains

    !> Calls the RKINTC, but makes sure that the input LABEL values do not get
    !! changed.
    subroutine rkintc_safe(l1, l2, l3, l4, l5, result)
        integer, intent(in), value :: l1, l2, l3, l4, l5
        real(real64), intent(out) :: result
        integer :: l1_local, l2_local, l3_local, l4_local, l5_local
        l1_local = l1
        l2_local = l2
        l3_local = l3
        l4_local = l4
        l5_local = l5
        call RKINTC(l1_local, l2_local, l3_local, l4_local, l5_local, result)
    end subroutine rkintc_safe

end function eval_matrixelement

subroutine eval_asfs(cimatrix)
    use grasp_kinds, only: real64, dp
    use orb_C
    use prnt_C
    use syma_C
    use eigv_C
    use eigvec1_C
    !use g2k_parameters
    !use matrixelement_common
    implicit none

    real(real64), dimension(:,:), intent(in) :: cimatrix

    integer :: k, i, j
    real(real64) :: expectation_value

    print *, "Hello World"
    print *, "NCF", NCF, "NW", NW
    !print *, NCFBL(1), NCFBL(2), NCFBL(3), NCFBL(4)
    print *, "NVEC", NVEC
    do k = 1, NVEC
        do i = 1, NCF
            print *, EVEC(i + (k - 1) * NCF), asf_coefficient(k, i)
        enddo
        print *
    enddo
    print *, "EVEC", EVEC(1), EVEC(2), EVEC(3)

    do k = 1, NVEC
        expectation_value = 0.0_dp
        do i = 1, NCF
            do j = 1, NCF
                expectation_value = expectation_value + asf_coefficient(k, i) * asf_coefficient(k, j) * cimatrix(i, j)
            enddo
        enddo
        print *, k, expectation_value
    enddo

contains

    !> Returns the i-th ASF coefficient of the k-th state
    function asf_coefficient(k, i)
        integer, intent(in) :: k, i
        real(real64) :: asf_coefficient
        asf_coefficient = EVEC(i + (k - 1) * NCF)
    end function asf_coefficient

end subroutine eval_asfs
