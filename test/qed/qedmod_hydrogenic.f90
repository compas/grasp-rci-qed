program qed_qedmod_hydrogenic_test
    use grasp_rciqed_kinds, only: real64, dp
    use parameter_def, only: NNN1, NNNP
    use orb_C, only: NW, NP, NAK
    use grasp_rciqed_qed_qedmod
    use grasp_rciqed_lib9290_init
    use grasptest_lib9290_setup, only: kappa_to_string
    implicit none

    type qedmod_values
        real(real64) :: overlap, uehling, wichkroll, selfenergy
    end type qedmod_values

    logical :: tests_passed = .true.
    integer :: n, k, l

    real(real64), parameter :: nuclear_z = 18.0_dp

    call lib9290_init_constants
    call lib9290_init_grid(nuclear_z)
    call init_nucleus ! Set up a Fermi nucleus
    call qedse_qedmod_init

    print *
    print *, "qedmod SE estimates for hydrogenic wavefunctions (DCBSRW routine vs SHAB_dirac)"
    print '(a12,"|",(a17,7(" "),"|"),4(a25,17(" "),"|"))', &
        '', 'Overlap', 'Orbital energy', 'qedmod SE', 'Uehling', 'Wichmann-Kroll'
    print '(a12," ",2a12," ",4(2a15,a12," "))', &
        '', 'Grasp', 'qedmod', &
        'Grasp', 'qedmod', 'diff', &
        'Grasp', 'qedmod', 'diff', &
        'Grasp', 'qedmod', 'diff', &
        'Grasp', 'qedmod', 'diff'

    do n = 1, 3
        do l = 0, n-1
            if(l>0) then
                call test_print_orbital(n, l)
            endif
            call test_print_orbital(n, -(l+1))
        end do
    end do

    if(.not.tests_passed) then
        print *, "qed_qedmod_hydrogenic_test: Tests failed."
        stop 1
    end if

contains

    subroutine init_nucleus
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: CVAC, C, PI, TENMAX, EXPMAX, EXPMIN, PRECIS, Z
        use npar_C, only: NPARM, PARM
        use nucpot_I

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        print *, TENMAX,EXPMAX,EXPMIN,PRECIS
        print *, CVAC, PI

        Z = nuclear_z

        NPARM = 2
        PARM(1) = 6.8839456274865651D-005
        PARM(2) = 9.8905913700962641D-006

        call NUCPOT
    end subroutine init_nucleus

    subroutine test_print_orbital(n, kappa)
        use grasptest_lib9290_setup, only: kappa_to_string
        use grasptest_testing

        integer, intent(in) :: n, kappa
        real(real64) :: energy_grasp, energy_qedmod
        type(qedmod_values) :: values_grasp, values_qedmod

        call qedse_hydrogenic_grasp(n, kappa, energy_grasp, values_grasp)
        call qedse_hydrogenic_qedmod(n, kappa, energy_qedmod, values_qedmod)

        print '(i5,a2,"(",i2,")  ",2es12.3," ",2es15.7,es12.3," ",2es15.7,es12.3," ",2es15.7,es12.3," ",2es15.7,es12.3)', &
            n, kappa_to_string(kappa), kappa, &
            values_grasp%overlap, values_grasp%overlap, &
            energy_grasp, energy_qedmod, energy_qedmod-energy_grasp, &
            values_grasp%selfenergy, values_qedmod%selfenergy, values_qedmod%selfenergy-values_grasp%selfenergy, &
            values_grasp%uehling, values_qedmod%uehling, values_qedmod%uehling-values_grasp%uehling, &
            values_grasp%wichkroll, values_qedmod%wichkroll, values_qedmod%wichkroll-values_grasp%wichkroll

        ! Note: tolerances are chosen so that the tests would pass.
        call test_isequal(tests_passed, "total energy", energy_grasp, energy_qedmod, 1e-5_dp)
        call test_isequal(tests_passed, "qedmod self-energy", values_grasp%selfenergy, values_qedmod%selfenergy, 1e-4_dp)
        call test_isequal(tests_passed, "Uehling VP", values_grasp%uehling, values_qedmod%uehling, 1e-4_dp)
        call test_isequal(tests_passed, "Wichmann-Kroll", values_grasp%wichkroll, values_qedmod%wichkroll, 1e-4_dp)
    end subroutine test_print_orbital

    subroutine qedse_hydrogenic_grasp(orbital_n, orbital_kappa, energy, values)
        use iso_fortran_env, only: error_unit
        use grasp_rciqed_kinds, only: real64, dp
        use parameter_def, only: NNN1, NNNP
        use def_C, only: Z
        use grid_C, only: R
        use dcbsrw_I
        use grasp_rciqed_qed_qedmod
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_cl=>cl, shab_z=>z, shab_r2=>r2, &
            shab_ii=>ii, shab_maxii=>maxii, shab_r=>r, shab_v=>v, shab_unuc=>unuc, &
            shab_knucl=>knucl, shab_rnucl=>rnucl, &
            shab_xr=>xr, shab_Gc=>Gc, shab_Fc=>Fc, shab_npoints=>npoints, &
            shab_nc=>nc, shab_KC=>Kc, shab_NnC=>NnC, &
            SHAB_wfunc_interpolate
        implicit none

        integer, intent(in) :: orbital_n, orbital_kappa
        real(real64), intent(out) :: energy
        type(qedmod_values), intent(out) :: values

        integer :: i
        real(real64) :: p(shab_maxii), q(shab_maxii), cp(shab_maxii), cq(shab_maxii)
        real(real64) :: ptemp, qtemp, selfenergy

        real(real64) :: orb_energy
        integer :: i_shab, i_g2k

        real(real64) :: SHAB_tint_args

        real(real64) :: dcwf_RG0, dcwf_RG(NNNP), dcwf_RF(NNNP)
        integer :: dcwf_MTP

        call DCBSRW(orbital_n, orbital_kappa, Z, energy, dcwf_RG0, dcwf_RG, dcwf_RF, dcwf_MTP)

        shab_nc = 1
        shab_KC(1) = orbital_kappa
        shab_NnC(1) = orbital_n
        shab_npoints(1) = dcwf_MTP
        do i = 1, dcwf_MTP
            shab_xr(i, 1) = R(i)
            shab_Gc(i, 1) = dcwf_RG(i)
            shab_Fc(i, 1) = dcwf_RF(i)
        end do

        do i = 1, shab_ii
            ! qedse_qedmod_init ensures that the qedmod grid lies within the same limits
            ! as the GRASP grid. This means that we can avoid extrapolation issues.
            call SHAB_wfunc_interpolate(orbital_kappa, orbital_n, shab_r(i), p(i), q(i))
            ! GRASP and QEDMOD have different conventions for Q(r). E.g. compare
            ! equation (15) in [QEDMOD, 2015] to the first equation in chapter 5.3
            ! of the [Grant, 2007] book.
            q(i) = -q(i)
        end do

        values = compute_qedmod_values(orbital_n, orbital_kappa, p, q)
    end subroutine qedse_hydrogenic_grasp

    subroutine qedse_hydrogenic_qedmod(n, kappa, energy, values)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: CVAC, PI, Z
        use grasp_rciqed_qed_qedmod, only: qedmod_initialized
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_cl=>cl, shab_z=>z, shab_r2=>r2, &
            shab_maxii=>maxii, shab_r=>r, shab_v=>v, shab_unuc=>unuc, &
            shab_knucl=>knucl, shab_rnucl=>rnucl, &
            shab_v_wk=>v_wk, shab_uehl=>uehl, &
            shab_tt=>tt, shab_aa=>aa, shab_cc=>cc, &
            shab_ii=>ii, &
            SHAB_dirac
        implicit none

        integer, intent(in) :: n, kappa
        real(real64), intent(out) :: energy
        type(qedmod_values), intent(out) :: values

        real(real64) :: p(shab_maxii), q(shab_maxii), cp(shab_maxii), cq(shab_maxii)
        real(real64) :: SHAB_sint, SHAB_tint, SHAB_tint_args
        real(real64) :: qedmodse_full

        if(.not.qedmod_initialized) then
            print *, "ERROR(libqed): qedmod shared data not initalized."
            print *, "qedse_qedmod_init() has to be called before the other routines can be used."
            stop 1
        end if

        ! Populate the p and q arrays with a hydrogenic wavefunction
        call SHAB_dirac(n, kappa, p, q)
        ! The energy of the hydrogenic wavefunction gets stored in p(ii+1), but is twice the value
        ! of the energy (?). Cf. SHAB_populate_hydrogenics()
        energy = 0.5_dp * p(shab_ii+1)
        values = compute_qedmod_values(n, kappa, p, q)
    end subroutine qedse_hydrogenic_qedmod

    function compute_qedmod_values(n, kappa, p, q) result(values)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: CVAC, PI, Z
        use grasp_rciqed_qed_qedmod, only: qedmod_initialized
        use grasp_rciqed_qed_qedmod_common, only: &
            shab_cl=>cl, shab_z=>z, shab_r2=>r2, &
            shab_maxii=>maxii, shab_r=>r, shab_v=>v, shab_unuc=>unuc, &
            shab_knucl=>knucl, shab_rnucl=>rnucl, &
            shab_v_wk=>v_wk, shab_uehl=>uehl, &
            shab_tt=>tt, shab_aa=>aa, shab_cc=>cc, &
            shab_ii=>ii, &
            SHAB_uehling, SHAB_tint, SHAB_sint, SHAB_wk, SHAB_se_pot_wav
        implicit none

        integer, intent(in) :: n, kappa
        real(real64), intent(in) :: p(shab_maxii), q(shab_maxii)
        type(qedmod_values) :: values

        real(real64) :: cp(shab_maxii), cq(shab_maxii)
        real(real64) :: magic_coef
        real(real64) :: SHAB_tint_args
        real(real64) :: qedmodse_full

        if(.not.qedmod_initialized) then
            print *, "ERROR(libqed): qedmod shared data not initalized."
            print *, "qedse_qedmod_init() has to be called before the other routines can be used."
            stop 1
        end if

        magic_coef = CVAC/PI*(Z/CVAC)**4/n**3

        ! Calculate the overlap
        values%overlap = SHAB_tint(0, p, q, p, q, shab_r, shab_v)

        ! Calculate Uehling potential VP contribution
        call SHAB_uehling(shab_maxii, shab_r, shab_unuc, shab_uehl)
        values%uehling = SHAB_sint(shab_uehl, p, q, p, q, shab_r, shab_v)

        ! Calculate Wichmann-Kroll VP contribution
        call SHAB_wk(shab_maxii, shab_r, shab_v_wk)
        values%wichkroll = SHAB_sint(shab_v_wk, p, q, p, q, shab_r, shab_v)*(CVAC**2)

        ! Calculate the self-energy
        call SHAB_se_pot_wav(n, kappa, shab_r, p, q, cp, cq)
        values%selfenergy = SHAB_tint(0, p, q, cp, cq, shab_r, shab_v)
    end function compute_qedmod_values

end program qed_qedmod_hydrogenic_test
