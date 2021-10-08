!> Private module containing common block definitions used in the QEDMOD code.
!!
!! A separate module is used so that it would be possible to use `implicit none`
!! in the main function body and also because the variable names clash with
!! GRASP global variables. `use ..., only: ... => ...` can be used to rename
!! the variables when using them to get around this.
!! The module __should not__ be used outside of the `se_qedmod.f90` file.
module grasp_rciqed_qed_qedmod_common
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit real*8 (a-h,o-z)

    include 'qedmod.inc'
    integer, parameter :: numRadPo = 10000, numwfMax = 10 ! have to agree with values in shabaev/wfunc_interpolate.f

    COMMON/SHAB_const/cl,z &
          /SHAB_grd/h,r1,r2,al,bt,ii &
          /SHAB_ncls/knucl,inucl,rnucl &
          /SHAB_pots/v(maxii),unuc(maxii),uehl(maxii),v_wk(maxii),  &
          vloc(maxii),vsemi_loc(maxii,6),wsemi_loc(maxii,6),ww(maxii) &
          /SHAB_ry/y(maxii),r(maxii) &
          /SHAB_fermi/aa,tt,cc,ro0,rnucl0 &
          /SHAB_DF/xr(numRadPo,numwfMax),Gc(numRadPo,numwfMax), &
          Fc(numRadPo,numwfMax),npoints(numwfMax),nc,KC(numwfMax),NnC(numwfMax)

    interface
        ! integr.f
        function SHAB_tint(k,p,q,a,b,r,v)
            include 'qedmod.inc'
            integer :: k
            real*8 p(maxii),q(maxii),a(maxii),b(maxii)
            real*8 r(maxii),v(maxii)
            real*8 SHAB_tint
        end function SHAB_tint

        function SHAB_sint(c,p,q,a,b,r,v)
            include 'qedmod.inc'
            real*8 c(maxii),p(maxii),q(maxii),a(maxii),b(maxii)
            real*8 r(maxii),v(maxii)
            real*8 SHAB_sint
        end function SHAB_sint

        ! uehling.f
        subroutine SHAB_uehling(maxii,r,unuc,uehl)
            integer :: maxii
            real*8 r(maxii),unuc(maxii),uehl(maxii)
        end subroutine SHAB_uehling

        ! dirac.f
        subroutine SHAB_dirac(n,kappa,p,q)
            include 'qedmod.inc'
            integer :: k, kappa
            real*8 p(maxii),q(maxii)
        end subroutine SHAB_dirac

        ! pot_wav.f
        subroutine SHAB_se_pot_wav(n,kappa,r,pp,qq,cp,cq)
            include 'qedmod.inc'
            integer :: k, kappa
            real*8 r(maxii),pp(maxii),qq(maxii),cp(maxii),cq(maxii)
        end subroutine SHAB_se_pot_wav

        ! wfunc_interpolate.f
        subroutine SHAB_wfunc_interpolate(KA,NnA,x,GA,FA)
            integer :: KA, NnA
            real*8 x, GA, FA
        end subroutine SHAB_wfunc_interpolate

        ! wk.f
        subroutine SHAB_wk(maxii,r,v_wk)
            integer :: maxii
            real*8 r(maxii),v_wk(maxii)
        end subroutine SHAB_wk
    end interface

end module grasp_rciqed_qed_qedmod_common

module grasp_rciqed_qed_qedmod
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

    logical :: qedmod_initialized = .false.

contains

    !> Estimates the QED self-energy of the `k`th orbital using the model operator
    !! approach described in [Shabev et al., 2013].
    !!
    !! The method relies on (modified) routines from the published implementation of
    !! the method [QEDMOD, 2015].
    !!
    !! ### References
    !!
    !!   - [Shabev et al., 2013] Shabaev, Tupitsyn and Yerokhin, Phys.Rev.A 88, 012513 (2013)
    !!   - [QEDMOD, 2015] Shabaev, Tupitsyn and Yerokhin, Comp.Phys.Comm 189, 175–181 (2015)
    !!
    function qedse_qedmod(k1, k2)
        use parameter_def, only: NNN1, NNNP
        use orb_C, only: NP, NH, NAK
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_maxii=>maxii, shab_r=>r, shab_v=>v
        implicit none

        integer, intent(in) :: k1, k2
        real(real64) :: qedse_qedmod

        integer :: i, orbital_n1, orbital_n2, orbital_kappa
        real(real64) :: p(shab_maxii), q(shab_maxii), cp(shab_maxii), cq(shab_maxii)

        real(real64) :: SHAB_tint

        if(.not.qedmod_initialized) then
            print *, "ERROR(libqed): Shabaev shared data not initalized."
            print *, "qedse_qedmod_init() has to be called before the other routines can be used."
            stop 1
        end if

        if(NAK(k1) /= NAK(k2)) then
            print '(a, ": (", i2, a2, " | ", i2, a2, ")")', &
                "WARNING(libqed:qedse_qedmod): Calculating SE matrix element for orbitals with differing kappas", &
                NP(k1), NH(k1), NP(k2), NH(k2)
            qedse_qedmod = 0.0_dp
            return
        endif
        orbital_kappa = NAK(k1)
        orbital_n1 = NP(k1)
        orbital_n2 = NP(k2)

        ! The QEDMOD code can only calculate the self-energy for s/p/d shells. So,
        ! for shells outside of this range, we define the self-energy contribution
        ! to be zero.
        !
        ! kappa for the f- shell is 3 and for the f shell is -4.
        if(orbital_kappa <= -4 .or. orbital_kappa >= 3) then
            qedse_qedmod = 0.0_dp
            return
        endif

        ! Calculate the self-energy. First, we'll interpolate the k1 orbital onto the
        ! p/q arrays. SHAB_se_pot_wav then populates the cp/cq arrays with the orbital
        ! that has been acted on by the self-energy operator. Then we'll interpolate
        ! the k2 orbital onto p/q and find the overlap of p/q and cp/cq. This yields
        ! the matrix element of the self-energy operator.
        call qedse_qedmod_wfinterpolate(k1, p, q)
        call SHAB_se_pot_wav(orbital_n1, orbital_kappa, shab_r, p, q, cp, cq)
        call qedse_qedmod_wfinterpolate(k2, p, q)
        qedse_qedmod = SHAB_tint(0, p, q, cp, cq, shab_r, shab_v)

    end function qedse_qedmod

    !> Interpolate the GRASP orbital `k` onto the QEDMOD grid.
    !!
    !! `p` and `q` are the interpolated output arrays on the QEDMOD grid.
    subroutine qedse_qedmod_wfinterpolate(k, p, q)
        use grid_C, only: N, R
        use npar_C, only: NPARM, PARM
        use orb_C, only: NP, NAK
        use wave_C, only: PF, QF
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_ii=>ii, shab_maxii=>maxii, shab_r=>r, shab_xr=>xr, &
            shab_Gc=>Gc, shab_Fc=>Fc, shab_npoints=>npoints, shab_nc=>nc, &
            shab_KC=>Kc, shab_NnC=>NnC
        implicit none

        integer, intent(in) :: k
        real(real64), intent(out) :: p(shab_maxii), q(shab_maxii)

        integer :: orbital_n, orbital_kappa
        integer :: i

        orbital_n = NP(k)
        orbital_kappa = NAK(k)

        ! SHAB_wfunc_interpolate(KA,NnA,x,GA,FA)
        !
        ! KA -- kappa, NnA -- n of the state we will interpolate
        ! x -- point to be interpolated
        ! GA, GA -- p, q values (output)
        !
        ! The function that is interpolated is read from SHAB_DF common block.
        !
        ! SHAB_DF
        !   /nc -- the number of orbitals stores
        !   /xr -- numRadPo array for RP values
        !   /Gc,Fc -- numRadPo,numwfMax arrays for PF/QF values
        !   /KC,NnC -- numwfMax arrays for kappa/n values of stored orbitals
        shab_nc = 1
        shab_KC(1) = orbital_kappa
        shab_NnC(1) = orbital_n
        shab_npoints(1) = N
        do i = 1, N
            shab_xr(i, 1) = R(i)
            shab_Gc(i, 1) = PF(i, k)
            shab_Fc(i, 1) = QF(i, k)
        end do

        do i = 1, shab_ii
            ! qedse_qedmod_init ensures that the Shabaev grid lies within the same limits
            ! as the Grasp2k grid. This means that we can avoid extrapolation issues.
            call SHAB_wfunc_interpolate(orbital_kappa, orbital_n, shab_r(i), p(i), q(i))
            ! G2K and QEDMOD have different conventions for Q(r). E.g. compare
            ! equation (15) in [QEDMOD, 2015] to the first equation in chapter 5.3
            ! of the [Grant, 2007] book.
            q(i) = -q(i)
        end do
    end subroutine qedse_qedmod_wfinterpolate

    subroutine qedse_qedmod_init
        use def_C, only: Z, CVAC, PI, FMTOAU
        use grid_C, only: N, R
        use npar_C, only: NPARM, PARM
        use grasp_rciqed_nucleus, only: fermi_rms
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_cl=>cl, shab_z=>z, shab_r2=>r2, &
            shab_maxii=>maxii, shab_r=>r, shab_v=>v, shab_unuc=>unuc, &
            shab_knucl=>knucl, shab_rnucl=>rnucl, &
            shab_v_wk=>v_wk, shab_uehl=>uehl, &
            shab_tt=>tt, shab_aa=>aa, shab_cc=>cc, &
            shab_ii=>ii
        implicit none

        if(qedmod_initialized) then
            print *, "WARNING: qedse_qedmod_init called multiple times."
            return
        endif

        ! Shabaev's routines write to file 11 -- we'll ignore that output
        open(unit=11, file='/dev/null', status='unknown')

        ! We need to set the speed of light (normally done in atom_data)
        shab_cl=CVAC

        ! The grid subroutine only needs r2 and z set.
        ! z is normally provided by the user, r2 is set by atom_data to 150/Z
        shab_z = Z
        shab_r2 = 150 / shab_z

        ! Set up the parameters for the Fermi nucleus.
        ! knucl determines the type of the nucleus: 0 -- point nucleus, 1 -- Fermi
        ! nucleus. In Grasp we look at the NPARM first to determine the type of the
        ! nuclear model.
        select case (NPARM)
        case (0)
            shab_knucl = 0
        case (2)
            shab_knucl = 1

            ! Stored in atomic units in both Grasp and Shabaev's code
            shab_cc = PARM(1)
            shab_aa = PARM(2)

            ! For consistency, we also set the t parameter, which we have to reconstruct
            ! from the a values. The two parameters are related by a = t/(4 log(3)).
            ! For some reason, Shabaev stores the t variable in femtometers (where it
            ! has the conventional value of 2.3 fm).
            shab_tt = shab_aa * 4.0_dp * log(3.0_dp) / FMTOAU

            ! PARM() stores the values in atomic units, and Shabaev's code assumes uses
            ! the same units for the input rnucl variable. ESTRMS appears to be agnostic
            ! towards the units of its inputs, as it should be.
            ! However, in Shabaev's code, rnucl = RMS * sqrt(5/3), hence the additional
            ! factor.
            shab_rnucl = fermi_rms(shab_aa, shab_cc) * sqrt(5.0_dp/3.0_dp)
        case default
            print *, "ERROR(libqed/qedse_qedmod_init): Invalid nuclear model."
            print '("NPAR/NPARM=",i0,", but only 0 or 2 allowed (for point or Fermi nucleus).")', NPARM
            stop 1
        end select

        ! Call the Shabaev initialization routines
        call SHAB_grid_args(1e-12_dp, R(N))
        call SHAB_nucl(shab_maxii, shab_r, shab_unuc, .false.)
        call SHAB_init_se

        ! Compute the hydrogenic wavefunctions, and both the local and non-local potential
        ! The arguments to SHAB_populate_hydrogenics control its verbosity.
        call SHAB_populate_hydrogenics(.false., 0)
        call SHAB_local_se_pot
        call SHAB_nonlocal_se_pot

        qedmod_initialized = .true.

    end subroutine qedse_qedmod_init

end module grasp_rciqed_qed_qedmod
