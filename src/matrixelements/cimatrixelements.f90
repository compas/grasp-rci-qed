!> Routines for evaluating the many-particle matrix elements of the Hamiltionian.
!!
!! Every routine normally has two methods: one that just takes the column and
!! row index of the CSF as an input, and another one which also accepts some
!! cached values (these variants are called `*_cached`).
!!
!! The cached versions are handy because normally the calculation happens in two
!! steps. In the first step, the angular parts are evaluated, and then in the
!! second step the one- and two-particle matrix elements are combined with the
!! angular information to form the complete matrix element. The first step,
!! however, is the same for all the operators of the same particle number and it
!! does not make sense to repeat it. The `*_cached` versions of the routines can
!! be used to achieve that.
!!
!! TODO: Actually, the two-particle elements (Coulomb and Breit) you can not split
!! up like that. They need their own `RKCO_GG` calls.
!!
!! TODO: Actually, coulomb and sms can share angular state.
module grasp_cimatrixelements
    use grasp_kinds, only: real64
    implicit none

    !> Calculates the many-particle CI matrix elements of the Dirac kinetic and
    !! nuclear potential parts of the Hamiltonian.
    interface dirac_potential
        module procedure dirac_potential, dirac_potential_cached
    end interface dirac_potential

    !> Matrix elements smaller than `cutoff` are set to zero.
    !!
    !! The use of this constant originates from the `matrix.f` file of the
    !! original `rci_mpi`.
    real(real64), parameter :: cutoff = 1d-20

contains

    ! TODO: single particle cache
    ! TODO: two particle cache

    !===========================================================================
    ! Routines for single-particle Dirac + nuclear potential matrix elements
    !---------------------------------------------------------------------------

    !> This also calls the `ONESCALAR` routine from lib92.
    function dirac_potential(ic, ir)
        use grasp_kinds, only: real64, dp
        use orb_C
        use onescalar_I
        implicit none

        integer, intent(in) :: ic, ir
        real(real64) :: dirac_potential

        real(real64) :: tshell(NW)
        integer :: k, l

        ! This part calculates the one-particle matrix elements. Based on
        ! setham_gg.f from rci_mpi.
        call ONESCALAR(ic, ir, k, l, tshell)
        dirac_potential = dirac_potential_cached(ic, ir, k, l, tshell)
    end function dirac_potential

    !> Calls the `IABINT` routine to evaluate the matrix element. The output from
    !! `ONESCALAR` can be passed as arguments.
    function dirac_potential_cached(ic, ir, k, l, tshell)
        use grasp_kinds, only: real64, dp
        use orb_C
        implicit none

        integer, intent(in) :: ic, ir, k, l
        real(real64), intent(in) :: tshell(NW)
        real(real64) :: dirac_potential_cached

        integer :: m
        real(real64) :: result

        ! It appears that when k==0, then there is no contribution? So, following
        ! setham_gg.f, we only do something when k /= 0.\
        dirac_potential_cached = 0.0_dp
        if(k /= 0 .and. k == l) then
            ! If k==l, then it seems that we have to add just the diagonal elements
            ! of the single particle matrix elements, weighted with TSHELL, to the
            ! matrix element.
            do m = 1, NW
                if(abs(tshell(m)) <= cutoff) cycle
                call iabint_safe(m, m, result)
                dirac_potential_cached = dirac_potential_cached + result * tshell(m)
            enddo
        elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
            call iabint_safe(k, l, result)
            dirac_potential_cached = dirac_potential_cached + result * tshell(1)
        endif
    end function dirac_potential_cached

    !> Calls the `IABINT` routine, but ensures that the input `k` and `l`
    !! arguments do not get swapped.
    subroutine iabint_safe(k, l, result)
        use grasp_kinds, only: real64
        use iabint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call IABINT(k, l, result)
    end

    !===========================================================================
    ! Routines for one-particle mass shift (normal mass shift; NMS) matrix
    ! elements
    !---------------------------------------------------------------------------

    function nms(ic, ir)
        use grasp_kinds, only: real64, dp
        use orb_C
        use onescalar_I
        implicit none

        integer, intent(in) :: ic, ir
        real(real64) :: nms

        real(real64) :: tshell(NW)
        integer :: k, l

        ! This part calculates the one-particle matrix elements. Based on
        ! setham_gg.f from rci_mpi.
        call ONESCALAR(ic, ir, k, l, tshell)
        nms = nms_cached(ic, ir, k, l, tshell)
    end function nms

    function nms_cached(ic, ir, k, l, tshell)
        use grasp_kinds, only: real64, dp
        use def_C, only: EMN
        use orb_C
        implicit none

        integer, intent(in) :: ic, ir, k, l
        real(real64), intent(in) :: tshell(NW)
        real(real64) :: nms_cached

        integer :: m
        real(real64) :: result, ATWINV

        ATWINV = 1.0_dp/EMN ! this is calculated in setham_gg.f

        ! The loops have the same logic as in the 1-particle DC.
        nms_cached = 0.0_dp
        if(k /= 0 .and. k == l) then
            do m = 1, NW
                if(abs(tshell(m)) <= cutoff) cycle
                call keint_safe(m, m, result)
                nms_cached = nms_cached + result * ATWINV * tshell(m)
            enddo
        elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
            call keint_safe(k, l, result)
            nms_cached = nms_cached + result * ATWINV * tshell(1)
        endif
    end function nms_cached

    !> Calls the `KEINT` routine, but ensures that the input `k` and `l`
    !! arguments do not get swapped.
    subroutine keint_safe(k, l, result)
        use grasp_kinds, only: real64
        use keint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call KEINT(k, l, result)
    end

    !===========================================================================
    ! Routines for vacuum polarization matrix elements
    !---------------------------------------------------------------------------

    function qed_vp(ic, ir)
        use grasp_kinds, only: real64, dp
        use orb_C
        use onescalar_I
        implicit none

        integer, intent(in) :: ic, ir
        real(real64) :: qed_vp

        real(real64) :: tshell(NW)
        integer :: k, l

        ! This part calculates the one-particle matrix elements. Based on
        ! setham_gg.f from rci_mpi.
        call ONESCALAR(ic, ir, k, l, tshell)
        qed_vp = qed_vp_cached(ic, ir, k, l, tshell)
    end function qed_vp

    function qed_vp_cached(ic, ir, k, l, tshell)
        use grasp_kinds, only: real64, dp
        use orb_C
        implicit none

        integer, intent(in) :: ic, ir, k, l
        real(real64), intent(in) :: tshell(NW)
        real(real64) :: qed_vp_cached

        integer :: m
        real(real64) :: result

        ! QED parts of the matrix element. The loops have the same logic as in
        ! the 1-particle DC.
        qed_vp_cached = 0.0_dp
        if(k /= 0 .and. k == l) then
            do m = 1, NW
                if(abs(tshell(m)) <= cutoff) cycle
                call vpint_safe(m, m, result)
                qed_vp_cached = qed_vp_cached + result * tshell(m)
            enddo
        elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
            call vpint_safe(k, l, result)
            qed_vp_cached = qed_vp_cached + result * tshell(1)
        endif
    end function qed_vp_cached

    !> Calls the `VPINT` routine, but ensures that the input `k` and `l`
    !! arguments do not get swapped.
    subroutine vpint_safe(k, l, result)
        use grasp_kinds, only: real64
        use vpint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call VPINT(k, l, result)
    end

    !===========================================================================
    ! Routines for Dirac + nuclear potential
    !---------------------------------------------------------------------------

    !===========================================================================
    ! Routines for two-particle Coulomb matrix elements
    !---------------------------------------------------------------------------

    function coulomb(ic, ir)
        use grasp_kinds, only: real64, dp
        use buffer_C, only: NVCOEF, COEFF, LABEL
        use cord_I
        use rkco_gg_I

        integer, intent(in) :: ic, ir
        real(real64) :: coulomb

        integer :: i
        real(real64) :: result

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
        ! The first `1` was called `incor`
        CALL RKCO_GG(ic, ir, CORD, 1, 1)

        ! Accumulate the matrix element
        coulomb = 0.0_dp
        do i = 1, NVCOEF
            if(abs(COEFF(i)) < cutoff) cycle
            call rkintc_safe(LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i), result)
            coulomb = coulomb + result * COEFF(i)
        enddo
    end function coulomb

    !===========================================================================
    ! Routines for two-particle mass shift (special mass shift; SMS)
    !---------------------------------------------------------------------------

    function sms(ic, ir)
        use grasp_kinds, only: real64, dp
        use def_C, only: EMN
        use buffer_C, only: NVCOEF, COEFF, LABEL
        use cord_I
        use rkco_gg_I
        use vint_I

        integer, intent(in) :: ic, ir
        real(real64) :: sms

        integer :: i
        real(real64) :: result1, result2, ATWINV

        ATWINV = 1.0_dp/EMN ! this is calculated in setham_gg.f

        NVCOEF = 0 ! this is done in SETHAM and then the variable used later
        ! The first `1` was called `incor`
        CALL RKCO_GG(ic, ir, CORD, 1, 1)

        sms = 0.0_dp
        do i = 1, NVCOEF
            if(abs(COEFF(i)) < cutoff) cycle
            if(LABEL(5,I) == 1) then
                call VINT(LABEL(1,i), LABEL(3,i), result1)
                call VINT(LABEL(2,i), LABEL(4,i), result2)
                sms = sms - result1 * result2 * ATWINV * COEFF(i)
            endif
        enddo
    end function sms

    !===========================================================================
    ! Routines for two-particle Breit matrix elements
    !---------------------------------------------------------------------------

    function breit(ic, ir)
        use grasp_kinds, only: real64, dp
        use buffer_C, only: NVCOEF, COEFF, LABEL
        use breid_I
        use rkco_gg_I
        use brint1_I
        use brint2_I
        use brint3_I
        use brint4_I
        use brint5_I
        use brint6_I

        integer, intent(in) :: ic, ir
        real(real64) :: breit

        integer :: itype, i
        real(real64) :: result

        ! The Breit interaction part of the matrix element
        breit = 0.0_dp
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

            if(itype > 0) then
                breit = breit + result * COEFF(i)
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
        enddo
    end function breit

    !===========================================================================
    ! Helper routines
    !---------------------------------------------------------------------------

    !> Calls the RKINTC, but makes sure that the input LABEL values do not get
    !! changed.
    subroutine rkintc_safe(l1, l2, l3, l4, l5, result)
        use grasp_kinds, only: real64
        use rkintc_I
        integer, value :: l1, l2, l3, l4, l5
        real(real64), intent(out) :: result
        call RKINTC(l1, l2, l3, l4, l5, result)
    end subroutine rkintc_safe

end module grasp_cimatrixelements
