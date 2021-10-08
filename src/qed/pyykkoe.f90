module grasp_rciqed_qed_pyykkoe
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

contains

    !> Calculates `(k1, k2)` matrix element of the effective local QED self-energy
    !! potential [Pyykkö & Zhao, 2003].
    !!
    !! Returns the value of the matrix element for the specified orbital pair `(k1, k2)`.
    !!
    !! It is assumed that the self-energy correction is given by a simple central scalar
    !! potential \f$ V(r) = B \exp(-\beta r^2) \f$. Pyykkö & Zhao then parametrize and fit
    !! the \f$B\f$ and \f$\beta\f$ values for \f$Z\f$-dependecy:
    !!
    !! \f[
    !!     B = b_0 + b_1 \cdot Z + b_2 \cdot Z^2,
    !!     \quad
    !!     \beta = \beta_0 + \beta_1 \cdot Z + \beta_2 \cdot Z^2
    !! \f]
    !!
    !! The values for the parameters are listed in Table 2. of the reference.
    !!
    !! The code is based on the implementation in GRASP92 by Christian Thierfelder,
    !! Peter Schwerdtfeger and Lukáš Félix Pašteka.
    !!
    !! ### References
    !!
    !!   - [Pyykkö & Zhao, 2003] Pyykkö and Zhao, J.Phys.B 36, 1469 (2003)
    !!
    function qedse_pyykkoe(k1, k2)
        use def_C, only: Z
        use grid_C, only: N, R, RP
        use tatb_C, only: MTP, TA
        use wave_C, only: PF, QF
        use quad_I
        implicit none

        integer :: k1, k2
        real(real64) :: qedse_pyykkoe

        integer :: i
        real(real64) :: B, beta

        ! Values for the parametrization of B and beta, from [Pyykkö & Zhao, 2003]
        real(real64), parameter :: b0 = -48.6166_dp
        real(real64), parameter :: b1 =  1.53666_dp
        real(real64), parameter :: b2 =  0.0301129_dp
        real(real64), parameter :: beta0 = -12751.3_dp
        real(real64), parameter :: beta1 =  916.038_dp
        real(real64), parameter :: beta2 =  5.7797_dp

        ! B and beta values for the current Z value
        B = b0 + b1*Z + b2*(Z**2)
        beta = beta0 + beta1*Z + beta2*(Z**2)

        ! We'll integrate f(r) = V(r) |ψ(r)|^2 , where ψ(r) is the Dirac wavefunction
        ! for the Pyykkö & Zhao potential V(r). The calculated expectation value gets
        ! returned.
        do i = 1, N
            TA(i) = (PF(i,k1)*PF(i,k2) + QF(i,k1)*QF(i,k2)) * B * exp(-beta*(R(i)**2)) * RP(i)
        end do
        MTP = N
        call QUAD(qedse_pyykkoe)
    end

end module grasp_rciqed_qed_pyykkoe
