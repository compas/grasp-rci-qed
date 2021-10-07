module grasp_rciqed_qed
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    implicit none

    integer, parameter :: nsetypes = 4
    character(*), dimension(nsetypes), parameter :: setypes_long = (/ &
        "Scaled hydrogenic (default)   ", &
        "QEDMOD (Shabaev et al., 2015) ", &
        "Flambaum & Ginges, 2005       ", &
        "Pyykkö & Zhao, 2003          "   & ! NOTE: ö is a 2-byte UTF-8 character
    /)
    character(*), dimension(nsetypes), parameter :: setypes_short = (/ &
        "hydrogenic", &
        "qedmod    ", &
        "flambaum  ", &
        "pyykkoe   "  &
    /)

contains

    subroutine init_vacuum_polarization
        use decide_C, only: LVP
        use grid_C, only: N, RP
        use ncdist_C, only: ZDIST
        use tatb_C, only: TB
        use vpilst_C, only: NVPI, FRSTVP
        use ncharg_I
        use vacpol_I
        use grasp_rciqed_qed_vp, only: vac2p4

        LVP = .TRUE.

        ! From AUXBLK. NCHARG set up the nuclear charge distrbution in ZDIST. If we're dealing
        ! with the PNC, then the routine actually does nothing, ZDIST stays zero and is
        ! ignored in VAC2 and VAC4.
        CALL NCHARG
        ! VACPOL calls VAC2 and VAC4, which populate TB with the (total) vacuum polarization
        ! potential.
        CALL VACPOL
        ! VACPOL puts the vacuum polarization potentials into TB, but the VPINT(F) routines
        ! use ZDSIT to actually evaluate the matrix elements.
        ZDIST(2:N) = TB(2:N)*RP(2:N)
        FRSTVP = .TRUE.
        NVPI = 0
        ! We'll also allocate the global array in grasp_rciqed_qed_vp and populate it with
        ! the potential value (without the RP multiplier, which is a QUAD-specific).
        allocate(vac2p4(N))
        vac2p4(1) = 0.0_dp
        vac2p4(2:N) = TB(2:N)
    end subroutine init_vacuum_polarization

    !> Populates the `matrix` with QED self-energy matrix elements for each orbital.
    !!
    !! `setype` determines the method used to estimate self-energy and can take the
    !! following values:
    !!
    !!   - `1` -- SE estimate based on hydrogenic wavefunctions (`qedse_mohr` in libqed)
    !!   - `2` -- Model operator approch for SE due to Shabaeva & Tupitsyna & Yerokhin (`qedse_shabaev` in libqed)
    !!   - `3` -- Flambaum & Ginges semi-empiric self-energy (`qedse_flambaum` in libqed)
    !!   - `4` -- Pyykkö & Zhao empiric self-energy (`qedse_pyykkoe` in libqed)
    !!
    !! `matrix` is assumed to be an `NW x NW` `real64` array.
    !!
    subroutine qedse(setype, matrix)
        use grasp_rciqed_kinds, only: real64, dp
        use orb_C, only: NW, NAK
        use grasp_rciqed_qed_pyykkoe
        use grasp_rciqed_qed_flambaum
        use grasp_rciqed_qed_qedmod

        integer :: setype
        real(real64) :: matrix(NW, NW)

        real(real64) :: seij, searray(NW)
        integer :: k, l
        real(real64) :: se_flam, se_flam_phi_l, se_flam_phi_f, se_flam_phi_g

        select case(setype)
        case(1)
            ! The original QED implementation, extrapolating from hydrogenic values,
            ! is still the default.
            call QED_SLFEN(searray)
            do k = 1, NW
                do l = k, NW
                    if(k == l) then
                        seij = searray(k)
                    else
                        seij = 0_dp
                    endif
                    matrix(k, l) = seij
                    matrix(l, k) = seij
                enddo
            enddo
        case(2)
            ! QEDMOD self-energy operator by Shabaev et al.
            call qedse_qedmod_init
            do k = 1, NW
                do l = k, NW
                    if(NAK(k) == NAK(l)) then
                        seij = qedse_qedmod(k, l)
                    else
                        seij = 0_dp
                    endif
                    matrix(k, l) = seij
                    matrix(l, k) = seij
                enddo
            enddo
        case(3)
            ! Self-energy operator by Flambaum & Ginges
            do k = 1, NW
                do l = k, NW
                    if(NAK(k) == NAK(l)) then
                        seij = qedse_flambaum(k, l, se_flam_phi_l, se_flam_phi_g, se_flam_phi_f)
                    else
                        seij = 0_dp
                    endif
                    matrix(k, l) = seij
                    matrix(l, k) = seij
                enddo
            enddo
        case(4)
            ! Self-energy operator by Pyykkö & Zhao
            do k = 1, NW
                do l = k, NW
                    if(NAK(k) == NAK(l)) then
                        seij = qedse_pyykkoe(k, l)
                    else
                        seij = 0_dp
                    endif
                    matrix(k, l) = seij
                    matrix(l, k) = seij
                enddo
            enddo
        case default
            stop "ERROR: Invalid setype for subroutine qedse"
        end select
    end subroutine qedse

end module grasp_rciqed_qed
