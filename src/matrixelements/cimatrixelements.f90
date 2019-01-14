!> Routines for evaluating the many-particle matrix elements of the Hamiltonian.
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
module grasp_rciqed_cimatrixelements
    use grasp_rciqed_kinds, only: real64, dp
    implicit none

    public dirac_potential, coulomb, coulomb_cached, breit, breit_split, &
        nms, sms, sms_cached, qed_vp, qed_se_mohr, qed_se
    public cutoff
    private

    !> Calculates the many-particle CI matrix elements of the Dirac kinetic and
    !! nuclear potential parts of the Hamiltonian.
    interface dirac_potential
        module procedure dirac_potential, dirac_potential_cached
    end interface dirac_potential

    !> Calculates the many-particle matrix element of the normal mass shift.
    interface nms
        module procedure nms, nms_cached
    end interface nms

    !> Calculates the many-particle matrix element of the QED vacuum polarization.
    interface qed_vp
        module procedure qed_vp, qed_vp_cached
    end interface qed_vp

    !> Calculates an estimate of the self-energy of a diagonal matrix element
    !! based on the Mohr hydrogenic values.
    interface qed_se_mohr
        module procedure qed_se_mohr, qed_se_mohr_cached
    end interface qed_se_mohr

    !> Calculates the self-energy many-body matrix element using the 1-particle
    !! matrix.
    interface qed_se
        module procedure qed_se, qed_se_cached
    end interface qed_se

    !> Matrix elements smaller than `cutoff` are set to zero.
    !!
    !! The use of this constant originates from the `setham_gg.f90` file.
    real(real64), parameter :: cutoff = 1e-12_dp

contains

    !===========================================================================
    ! Routines for single-particle Dirac + nuclear potential matrix elements
    !---------------------------------------------------------------------------

    !> This also calls the `ONESCALAR` routine from lib9290.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The Dirac kinetic + nuclear potential expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{D}} + \hat{V}_{\textrm{nucl.}}|\Psi_{\textrm{ir}}\rangle\f$.
    function dirac_potential(ic, ir)
        use grasp_rciqed_kinds, only: real64, dp
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
    !! `ONESCALAR` must be passed as arguments.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @param k, l, tshell The output from `ONESCALAR`.
    !! @returns The Dirac kinetic + nuclear potential expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{D}} + \hat{V}_{\textrm{nucl.}}|\Psi_{\textrm{ir}}\rangle\f$.
    function dirac_potential_cached(ic, ir, k, l, tshell)
        use grasp_rciqed_kinds, only: real64, dp
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
        use grasp_rciqed_kinds, only: real64
        use iabint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call IABINT(k, l, result)
    end

    !===========================================================================
    ! Routines for one-particle mass shift (normal mass shift; NMS) matrix
    ! elements
    !---------------------------------------------------------------------------

    !> This also calls the `ONESCALAR` routing from lib9290.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The normal mass shift expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{NMS}}|\Psi_{\textrm{ir}}\rangle\f$.
    function nms(ic, ir)
        use grasp_rciqed_kinds, only: real64, dp
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

    !> The `ONESCALAR` output has to be passed as arguments here.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @param k, l, tshell The output from `ONESCALAR`.
    !! @returns The normal mass shift expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{NMS}}|\Psi_{\textrm{ir}}\rangle\f$.
    function nms_cached(ic, ir, k, l, tshell)
        use grasp_rciqed_kinds, only: real64, dp
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
        use grasp_rciqed_kinds, only: real64
        use keint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call KEINT(k, l, result)
    end

    !===========================================================================
    ! Routines for vacuum polarization matrix elements
    !---------------------------------------------------------------------------

    !> This also calls the `ONESCALAR` routing from lib9290.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The QED vacuum polarization expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{VP}}|\Psi_{\textrm{ir}}\rangle\f$.
    function qed_vp(ic, ir)
        use grasp_rciqed_kinds, only: real64, dp
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

    !> The `ONESCALAR` output has to be passed as arguments here.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @param k, l, tshell The output from `ONESCALAR`.
    !! @returns The QED vacuum polarization expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{VP}}|\Psi_{\textrm{ir}}\rangle\f$.
    function qed_vp_cached(ic, ir, k, l, tshell)
        use grasp_rciqed_kinds, only: real64, dp
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
        use grasp_rciqed_kinds, only: real64
        use vpint_I
        integer, value :: k, l
        real(real64), intent(out) :: result
        call VPINT(k, l, result)
    end

    !===========================================================================
    ! Routines for two-particle Coulomb matrix elements
    !---------------------------------------------------------------------------

    !> Also calls the `RKCO_GG` routine.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The Coulomb interaction expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{C}}|\Psi_{\textrm{ir}}\rangle\f$.
    function coulomb(ic, ir)
        use grasp_rciqed_kinds, only: real64
        use buffer_C, only: NVCOEF
        use cord_I
        use rkco_gg_I

        integer, intent(in) :: ic, ir
        real(real64) :: coulomb

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
        coulomb = coulomb_cached(ic, ir)
    end function coulomb

    !> Assumes that `buffer_C` is set up properly by `RKCO_GG`.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The Coulomb interaction expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{C}}|\Psi_{\textrm{ir}}\rangle\f$.
    function coulomb_cached(ic, ir)
        use grasp_rciqed_kinds, only: real64
        use buffer_C, only: NVCOEF, COEFF, LABEL

        integer, intent(in) :: ic, ir
        real(real64) :: coulomb_cached

        integer :: i
        real(real64) :: result

        ! Accumulate the matrix element
        coulomb_cached = 0.0_dp
        do i = 1, NVCOEF
            if(abs(COEFF(i)) < cutoff) cycle
            call rkintc_safe(LABEL(1,i), LABEL(2,i), LABEL(3,i), LABEL(4,i), LABEL(5,i), result)
            coulomb_cached = coulomb_cached + result * COEFF(i)
        enddo
    end function coulomb_cached

    !> Calls the RKINTC, but makes sure that the input LABEL values do not get
    !! changed.
    subroutine rkintc_safe(l1, l2, l3, l4, l5, result)
        use grasp_rciqed_kinds, only: real64
        use rkintc_I
        integer, value :: l1, l2, l3, l4, l5
        real(real64), intent(out) :: result
        call RKINTC(l1, l2, l3, l4, l5, result)
    end subroutine rkintc_safe

    !===========================================================================
    ! Routines for two-particle mass shift (special mass shift; SMS)
    !---------------------------------------------------------------------------

    !> Also calls the `RKCO_GG` routine.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The special mass shift expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{SMS}}|\Psi_{\textrm{ir}}\rangle\f$.
    function sms(ic, ir)
        use grasp_rciqed_kinds, only: real64
        use buffer_C, only: NVCOEF
        use cord_I
        use rkco_gg_I

        integer, intent(in) :: ic, ir
        real(real64) :: sms

        NVCOEF = 0 ! this is done in SETHAM and then the variable used later
        ! The first `1` was called `incor`
        CALL RKCO_GG(ic, ir, CORD, 1, 1)
        sms = sms_cached(ic, ir)
    end function sms

    !> Assumes that `buffer_C` is set up properly by `RKCO_GG`.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The special mass shift expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{SMS}}|\Psi_{\textrm{ir}}\rangle\f$.
    function sms_cached(ic, ir)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: EMN
        use buffer_C, only: NVCOEF, COEFF, LABEL
        use vint_I

        integer, intent(in) :: ic, ir
        real(real64) :: sms_cached

        integer :: i
        real(real64) :: result1, result2, ATWINV

        ATWINV = 1.0_dp/EMN ! this is calculated in setham_gg.f

        sms_cached = 0.0_dp
        do i = 1, NVCOEF
            if(abs(COEFF(i)) < cutoff) cycle
            if(LABEL(5,I) == 1) then
                call VINT(LABEL(1,i), LABEL(3,i), result1)
                call VINT(LABEL(2,i), LABEL(4,i), result2)
                sms_cached = sms_cached - result1 * result2 * ATWINV * COEFF(i)
            endif
        enddo
    end function sms_cached

    !===========================================================================
    ! Routines for two-particle Breit matrix elements
    !---------------------------------------------------------------------------

    !> Returns the matrix element of the Breit operator.
    function breit(ic, ir)
        use grasp_rciqed_kinds, only: real64, dp

        integer, intent(in) :: ic, ir
        real(real64) :: breit

        real(real64) :: breit_core, breit_noncore, elsto

        if(ic == 1 .and. ir == 1) then
            call breit_split(ic, ir, breit_core, breit_noncore)
            breit = breit_core + breit_noncore
        elseif(ic == ir) then
            call breit_split( 1,  1, elsto, breit_noncore)
            call breit_split(ic, ir, breit_core, breit_noncore)
            if(breit_core /= 0.0_dp) then
                ! NOTE: comparing floats here with equality, but `breit_split`
                ! should set it to _exactly_ 0.0_dp and it should not be modified
                ! at all, unless the value is significant, so this should work
                ! as expected.
                print '("┌ WARNING: Unexpected non-zero value for breit_core")'
                print '("│ ir, ic = ",i0,", ", i0)', ir, ic
                print '("└ @ grasp_rciqed_cimatrixelements % breit W101")'
            endif
            breit = elsto + breit_noncore
        else ! ic /= ir
            call breit_split(ic, ir, breit_core, breit_noncore)
            if(breit_core /= 0.0_dp) then
                ! NOTE: comparing floats here with equality, but `breit_split`
                ! should set it to _exactly_ 0.0_dp and it should not be modified
                ! at all, unless the value is significant, so this should work
                ! as expected.
                print '("┌ WARNING: Unexpected non-zero value for breit_core")'
                print '("│ ir, ic = ",i0,", ", i0)', ir, ic
                print '("└ @ grasp_rciqed_cimatrixelements % breit W101")'
            endif
            breit = breit_noncore
        endif
    end function breit

    !> Calculates the Breit matrix element, but calculates the core part only
    !! sometimes.
    !!
    !! The `breit_core` argument is set to non-zero only if `ic == ir == 1`.
    !!
    !! This behaviour was reverse-engineered from `setham_gg.f90`. It appears
    !! that the contribution to the diagnoal elements from the core electrons
    !! is only calculated `ic == ir == 1`. For those contributions, `LABEL(6,i)`
    !! is set to a negative value and in `setham_gg` they are accumulated into
    !! `ELSTO`, not to `ELEMNT`.
    !!
    !! The assumption is that this contribution is the same for all the diagonal
    !! elements and therefore does not have to be calculated again, probably to
    !! optimize the calculation of the Breit matrix elements.
    !!
    !! This behaviour is completely undocumented, however, and the API does not
    !! give any hints. It just appears that `RKCI_GG` / `BREID` populate the
    !! `LABEL` arrays differently for some `ic`/`ir` values.
    subroutine breit_split(ic, ir, breit_core, breit_noncore)
        use grasp_rciqed_kinds, only: real64, dp
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
        real(real64), intent(out) :: breit_core, breit_noncore

        integer :: itype, i
        real(real64) :: result

        ! The Breit interaction part of the matrix element
        breit_core = 0.0_dp
        breit_noncore = 0.0_dp
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

            ! Note: LABEL(6,i) can also be negative. See breid.f90.
            !
            ! In breid.f90, the sign of ITYPE is controlled by ISG, which is
            ! determined by:
            !
            !     ISG = 1
            !     IF (JA == JB) THEN
            !        IF (ICORE(IA1)/=0 .AND. ICORE(IB1)/=0) THEN
            !           IF (JA > 1) RETURN
            !           ISG = -1
            !        ENDIF
            !     ENDIF
            !
            ! This means that it can sometimes also be negative. It seems to be
            ! negative if the value is for a contributions from the filled core
            ! shells. It would make sense that these contributions are the same
            ! for all CSFs, and hence can be "shared" for the diagonal elements.
            !
            ! This is what appears to be happening with ELSTO -- it stores the
            ! common part of the diagonal for each block.
            !
            ! From the original setham_gg routine:
            !
            !    IF (LABEL(6,I) > 0) THEN
            !       ELEMNT = ELEMNT + CONTR
            !    ELSE
            !    !            ...It comes here only when ic=ir=1
            !    !               clue: rkco<-breid<-talk<-label(6,i)
            !       NCORE = NCORE + 1
            !       ELSTO = ELSTO + CONTR
            !       write(fh_hmat, '(" % setting ELSTO")')
            !    ENDIF
            !
            ! So, to correctly calculate the matrix element, we should also
            ! add the ones that normally go to ELSTO to the Breit matrix element.
            if(LABEL(6,i) > 0) then
                breit_noncore = breit_noncore + result * COEFF(i)
            elseif(LABEL(6,i) < 0) then
                if(ir /= 1 .or. ir /= 1) then
                    ! The assumption is that `ELSTO` is only set when ic == ir == 1.
                    ! To make sure we catch the error if that assumption turns out
                    ! to be false, we will just crash the program here.
                    print '("┌ ERROR: Unexpected negative LABEL(6,i)")'
                    print '("│ ir, ic = ",i0,", ", i0)', ir, ic
                    print '("└ @ grasp_rciqed_cimatrixelements % breit_split E201")'
                    error stop
                endif
                breit_core = breit_core + result * COEFF(i)
            else
                ! Similar to above -- the assumption here is that `LABEL(6,:)`
                ! will never be zero.
                print '("┌ ERROR: Unexpected zero LABEL(6,i)")'
                print '("│ ir, ic = ",i0,", ", i0)', ir, ic
                print '("└ @ grasp_rciqed_cimatrixelements % breit_split E202")'
                error stop
            endif
        enddo
    end subroutine breit_split

    !===========================================================================
    ! Routines for calculating self-energy contributions for _diagonal_ elements
    ! based on hydrogenic SE values tabulated by Mohr et.al.
    !---------------------------------------------------------------------------

    !> This version initializes the `slfint` array, which contains the per-orbital
    !! self-energy contributions.
    !!
    !! @param ic The index of the CSF.
    !! @returns The self-energy expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{SE}}|\Psi_{\textrm{ic}}\rangle\f$.
    function qed_se_mohr(ic)
        use grasp_rciqed_kinds, only: real64, dp
        use parameter_def, only: NNNW
        use qed_slfen_I
        implicit none

        integer, intent(in) :: ic
        real(real64) :: qed_se_mohr

        real(real64) :: slfint(NNNW)

        ! Based on matrix.f90
        call QED_SLFEN(slfint)
        qed_se_mohr = qed_se_mohr_cached(ic, slfint)
    end function qed_se_mohr

    !> Calculates the QED contribution to a CI diagonal element using the orbital
    !! self-energy values in the `slfint` array.
    !!
    !! @param ic The index of the CSF.
    !! @param slfint A `NW`-element array containing the per-orbital self-energy
    !!   values.
    !! @returns The self-energy expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{SE}}|\Psi_{\textrm{ic}}\rangle\f$.
    function qed_se_mohr_cached(ic, slfint)
        use grasp_rciqed_kinds, only: real64, dp
        use parameter_def, only: NNNW
        use orb_C
        use iq_I
        implicit none

        integer, intent(in) :: ic
        real(real64), intent(in) :: slfint(NNNW)
        real(real64) :: qed_se_mohr_cached

        integer :: i

        ! Based on matrix.f90
        qed_se_mohr_cached = 0.0D00
        do i = 1, NW
            qed_se_mohr_cached = qed_se_mohr_cached + IQ(I,IC) * slfint(I)
        end do
    end function qed_se_mohr_cached

    !===========================================================================
    ! Generic routines for self-energy, based on the 1-particle matrix.
    !---------------------------------------------------------------------------

    !> This also calls the `ONESCALAR` routing from lib9290.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @returns The QED self-energy expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{VP}}|\Psi_{\textrm{ir}}\rangle\f$.
    function qed_se(sematrix, ic, ir)
        use grasp_rciqed_kinds, only: real64, dp
        use orb_C
        use onescalar_I
        implicit none

        real(real64), intent(in) :: sematrix(NW, NW)
        integer, intent(in) :: ic, ir
        real(real64) :: qed_se

        real(real64) :: tshell(NW)
        integer :: k, l

        ! This part calculates the one-particle matrix elements. Based on
        ! setham_gg.f from rci_mpi.
        call ONESCALAR(ic, ir, k, l, tshell)
        qed_se = qed_se_cached(sematrix, ic, ir, k, l, tshell)
    end function qed_se

    !> The `ONESCALAR` output has to be passed as arguments here.
    !!
    !! @param ic, ir The row and column indices of the CSFs.
    !! @param k, l, tshell The output from `ONESCALAR`.
    !! @returns The QED self-energy expectation value
    !!   \f$\langle\Psi_{\textrm{ic}}|\hat{H}_{\textrm{VP}}|\Psi_{\textrm{ir}}\rangle\f$.
    function qed_se_cached(sematrix, ic, ir, k, l, tshell)
        use grasp_rciqed_kinds, only: real64, dp
        use orb_C
        implicit none

        real(real64), intent(in) :: sematrix(NW, NW)
        integer, intent(in) :: ic, ir, k, l
        real(real64), intent(in) :: tshell(NW)
        real(real64) :: qed_se_cached

        integer :: m
        real(real64) :: result

        ! QED parts of the matrix element. The loops have the same logic as in
        ! the 1-particle DC.
        qed_se_cached = 0.0_dp
        if(k /= 0 .and. k == l) then
            do m = 1, NW
                if(abs(tshell(m)) <= cutoff) cycle
                qed_se_cached = qed_se_cached + sematrix(m, m) * tshell(m)
            enddo
        elseif(k /= 0 .and. k /= l .and. abs(tshell(1)) > cutoff) then
            qed_se_cached = qed_se_cached + sematrix(k, l) * tshell(1)
        endif
    end function qed_se_cached

end module grasp_rciqed_cimatrixelements
