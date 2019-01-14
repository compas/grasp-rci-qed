module grasp_rciqed_qed_flambaum
    implicit none

contains

    !> Estimates the QED self-energy of the `k`th orbital using potentials presented
    !! in [Flambaum & Ginges, 2005].
    !!
    !! Returns the self-energy estimate for the orbital `k`. Additional output will
    !! stored in the `phi_l`, `phi_f`, `phi_g` variables, which will contain the
    !! corresponding parts of the self-energy estimate.
    !!
    !! In the article the self-energy is assumed to be a spherically symmetric scalar
    !! potential, which has further been split up into three terms must be summed up. I.e.
    !! we're calculating equation (11), but without the Uehling a WC bits:
    !!
    !! \f[
    !!     \Phi_{\textrm{rad}}(r) = \Phi_l(r) + \Phi_f(r) + \Phi_g(r)
    !! \f]
    !!
    !! The code is based on the implementation in GRASP92 by Christian Thierfelder,
    !! Peter Schwerdtfeger and Lukáš Félix Pašteka.
    !!
    !! ### References
    !!
    !! [Flambaum & Ginges, 2005] Flambaum and Ginges, Phys.Rev.A 72, 052115 (2005)
    !! [Thierfelder & Schwerdtfeger, 2010] Thierfelder and Schwerdtfeger, Phys.Rev.A 82, 062503 (2010)
    !!
    function qedse_flambaum(k1, k2, phi_l, phi_f, phi_g)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: Z, CVAC, PI
        use grid_C, only: N, R, RP
        use npar_C, only: NPARM, PARM
        use npot_C, only: ZZ
        use orb_C, only: NP
        use tatb_C, only: MTP, TA
        use wave_C, only: PF, QF
        use quad_I

        integer, intent(in) :: k1, k2
        real(real64) :: qedse_flambaum
        real(real64), intent(out) :: phi_l, phi_f, phi_g

        integer :: i, NT
        real(real64) :: result, BZ, AZ, afit, g1, g2, diffnucpot, x, xi

        ! We need to perform inner integrations for Phi_g and Phi_f, which are done numerically.
        ! We do this on the following single-variable grid of NTMAX points. tt will contain the
        ! function/integrand values and ttx are the gridpoints.
        integer, parameter :: NTMAX = 10000
        real(real64) :: tt(NTMAX), ttx(NTMAX)

        ! Phi_l contribution, eq (9)
        ! -------------------------------------------------------------------------
        BZ = 0.074_dp + 0.35_dp*Z/CVAC
        do i = 1, N
            TA(i) = (PF(i,k1)*PF(i,k2) + QF(i,k1)*QF(i,k2)) * exp(-Z*R(i)) * RP(i)
        end do
        MTP = N
        call QUAD(result)
        phi_l = BZ*(Z/CVAC)*(Z/CVAC)*(Z/CVAC)*Z * result


        ! Phi_f contribution, eq (10)
        ! -------------------------------------------------------------------------

        ! Calculating the polynomial part of A(Z, r) and storing it in afit.

        if(k1 == k2) then
            ! Special fit for the 1s and 2s orbitals, recommended by Flambaum. But this
            ! can only be used for the diagonal elements. The fitting etc. is described
            ! in [Thierfelder & Schwerdtfeger, 2010].
            if      (NP(k1) == 1) then
                afit = 0.7645_dp+0.00230_dp*Z/(1+exp((Z/112.930_dp)**5))
            else if (NP(k1) == 2) then
                afit = 0.7912_dp+0.00629_dp*Z/(1+exp((Z/101.636_dp)**5))
            else if (NP(k1) == 3) then
                afit = 0.7980_dp+0.00738_dp*Z/(1+exp((Z/101.611_dp)**5))
            else if (NP(k1) == 4) then
                afit = 0.8009_dp+0.00779_dp*Z/(1+exp((Z/101.047_dp)**5))
            else if (NP(k1) == 5) then
                afit = 0.8023_dp+0.00799_dp*Z/(1+exp((Z/101.632_dp)**5))
            else if (NP(k1) == 6) then
                afit = 0.8032_dp+0.00813_dp*Z/(1+exp((Z/101.607_dp)**5))
            else
                afit = 0.8037_dp+0.00824_dp*Z/(1+exp((Z/101.591_dp)**5))
            end if
        else
            x = (Z-80)/CVAC
            afit = ((0.169*x - 2.128)*x - 1.976)*x*x + 1.071
        endif

        do i = 1, N
            AZ = afit * R(i) * CVAC / (R(i) * CVAC + 0.07_dp * (Z/CVAC)**2)

            ! t-integration
            call se_t_grid(tt, ttx, NTMAX, NT, R(i), 'F ')
            call quad_t_us2(tt, ttx, NTMAX, NT, RESULT)

            ! ZZ(i)/R(i) gives the nuclear potential, because ZZ(i) is the (effective)
            ! nuclear charge
            TA(i) = (PF(i,k1)*PF(i,k2) + QF(i,k1)*QF(i,k2)) * RESULT * (ZZ(i) / R(i)) * AZ * RP(i)
        end do
        MTP = N
        call QUAD(result)
        phi_f = result/(CVAC*PI)


        ! Phi_g contribution, eq (7)
        ! -------------------------------------------------------------------------
        do i = 1, N
            ! t-integration -- we need to do it twice because of the derivative
            ! Once the expression in the article and then also its derivative.
            call se_t_grid(tt, ttx, NTMAX, NT, R(i), 'G1')
            call quad_t_us2(tt, ttx, NTMAX, NT, result)

            ! Finding the derivative of the nuclear potential
            if(NPARM == 0) then
                ! point nucleus
                diffnucpot = ZZ(i)/R(i)**2
            else
                ! Gaussian nucleus
                xi = 1.5_dp / PARM(1)**2
                diffnucpot = ZZ(I)/R(i)**2 - 2*Z*sqrt(xi/PI)/R(i)*exp(-xi*R(i)**2)
            end if

            g1 = result*diffnucpot

            ! The second t-integration
            call se_t_grid(tt, ttx, NTMAX, NT, R(i), 'G2')
            call quad_t_us2(tt, ttx, NTMAX, NT, result)
            ! TODO: There might be a potential sign error on the following line (where did the
            ! minus from the exponential go to that should've come downstairs with the derivative
            ! as well?)
            g2 = 2 * result * CVAC * ZZ(i) / R(i)

            TA(i) = (PF(i,k1)*QF(i,k2) + QF(i,k1)*PF(i,k2)) * (g1 + g2 - diffnucpot) * RP(i)
        end do
        MTP = N
        call QUAD(result)
        phi_g = result / (4 * CVAC**2 * PI)

        qedse_flambaum = phi_l + phi_f + phi_g
    end function qedse_flambaum

    !> Calculate the integral-kernel for the radiative \f$\Phi_f\f$ potential.
    !!
    !! Equation 10 in [Flambaum & Ginges, 2005].
    !!
    function phi_f_inner(t, r)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: Z, CVAC, PI

        real(real64) :: t, r, phi_f_inner
        real(real64) :: x1, x2, x3, x4, x6, x7

        x1 = 1/sqrt(t**2 - 1)
        x2 = 1 - 1/(2 * t**2)
        x3 = log(t**2 - 1)
        x4 = 4*log(CVAC/Z + 0.5_dp)
        x6 = 1/t**2
        x7 = exp(-2*t*r*CVAC)
        phi_f_inner = x1*(x2*(x3 + x4) + x6 - 1.5_dp)*x7
    end function phi_f_inner

    !> Calculate the first part of the integral-kernel for the radiative \f$\Phi_g\f$
    !! potential.
    !!
    !! This is the first term from the product rule in equation (7) of [Flambaum & Ginges, 2005].
    !!
    function phi_g1_inner(t,r)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: CVAC, PI

        real(real64) :: t, r, phi_g1_inner

        phi_g1_inner = exp(-2*t*r*CVAC) / ((t**2) * sqrt(t**2 - 1))
    end function phi_g1_inner

    !> Calculate the second part of the integral-kernel for the radiative \f$\Phi_g\f$
    !! potential.
    !!
    !! This is the second term from the product rule in equation (7) of [Flambaum & Ginges, 2005].
    !!
    function phi_g2_inner(t,r)
        use grasp_rciqed_kinds, only: real64, dp
        use def_C, only: CVAC, PI

        real(real64) :: t, r, phi_g2_inner

        phi_g2_inner = exp(-2*t*r*CVAC)/(t * sqrt(t**2 - 1))
    end function phi_g2_inner

    !> Uses the Simpson's rule for a non-uniform grid `t` to evaluate the integral of `f`.
    !!
    !! `t` must be an array of length `ntmax` and will contain the gridpoints. `f` is an
    !! array of the same lengt, containing the values of the integrand for each `t(i)`.
    !! The integral is then evaluated between `t(1)` and `t(nt)` and stored in `result`.
    !!
    subroutine quad_t_us2(f, t, ntmax, nt, result)
        use grasp_rciqed_kinds, only: real64, dp

        integer :: ntmax, nt
        real(real64) :: f(ntmax), t(ntmax), result

        integer :: J1, J2, J3, j
        real(real64) :: f1, f2, f3, x1, x2, x3, t1, t2, t3, sum

        ! Find the starting point for the integration. This is the point where the
        ! integrand TA first becomes significant, and is chosen to prevent division
        ! by zero.
        sum = 0.0_dp
        do j = 1, nt-3, 2
            f1=f(j); f2=f(j+1); f3=f(j+2)
            t1=t(j); t2=t(j+1); t3=t(j+2)

            x1 = f1*(t2-t3)*(2*t1-3*t2+t3)
            x2 = f2*(t1-t3)**2
            x3 = f3*(t1-t2)*(t1-3*t2+2*t3)

            sum = sum - (t1-t3)*(x1+x2-x3)/(6*(t1-t2)*(t2-t3))
        end do
        result = sum
    end

    !> Populates the `t` and `f` arrays (of length `ntmax`) with grid (`t`) and integrand
    !! (`f`) values to prepare for the inner integrations of the radiative potentials.
    !!
    !! * `nt` will contain the largest array index which still contain non-neglible
    !!   integrand values.
    !! * `r` is the radius for which the integral kernels will be integrated.
    !! * `nterm` is is a string that determines which integrand will be used to populate
    !!   `f` (possible values are `F`, `G1` or `G2`).
    !!
    subroutine se_t_grid(f, t, ntmax, nt, r, nterm)
        use grasp_rciqed_kinds, only: real64, dp

        integer :: j, ntmax, nt
        character*2 :: nterm
        real(real64) :: f(ntmax), t(ntmax), r

        do j = 1, ntmax
            t(j) = 1.0_dp + 10.0_dp**(-(12.0_dp - j * 0.01_dp))

            if(nterm == 'F ') then
                f(j) = phi_f_inner(t(j), r)
            else if(nterm == 'G1') then
                f(j) = phi_g1_inner(t(j), r)
            else if(nterm == 'G2') then
                f(j) = phi_g2_inner(t(j), r)
            else
                print *, "ERROR: Bad nterm passed to se_t_grid()"
                print *, "nterm = '", nterm, "', expected 'F ', 'G1' or 'G2'"
                print *, "in rci/se_flam.90"
                stop 1
            end if

            if(abs(f(j)) <= 1.0e-32_dp) then
                nt = j
                exit
            end if
        end do
    end subroutine se_t_grid

end module grasp_rciqed_qed_flambaum
