module grasp_rciqed_qed
    implicit none

contains

    subroutine init_vacuum_polarization
        use decide_C, only: LVP
        use grid_C, only: N, RP
        use ncdist_C, only: ZDIST
        use tatb_C, only: TB
        use vpilst_C, only: NVPI, FRSTVP
        use ncharg_I
        use vacpol_I

        LVP = .TRUE.

        ! From AUXBLK
        CALL NCHARG
        CALL VACPOL
        ZDIST(2:N) = TB(2:N)*RP(2:N)
        FRSTVP = .TRUE.
        NVPI = 0
    end subroutine init_vacuum_polarization

    !> Populates the `matrix` with QED self-energy matrix elements for each orbital.
    !!
    !! `setype` determines the method used to estimate self-energy and can take the
    !! following values:
    !!
    !!   - `0` -- SE estimate based on hydrogenic wavefunctions (`qedse_mohr` in libqed)
    !!   - `1` -- Model operator approch for SE due to Shabaeva & Tupitsyna & Yerokhin (`qedse_shabaev` in libqed)
    !!   - `2` -- Flambaum & Ginges semi-empiric self-energy (`qedse_flambaum` in libqed)
    !!   - `3` -- Pyykkö & Zhao empiric self-energy (`qedse_pyykkoe` in libqed)
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
        case(0)
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
        case(1)
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
        case(2)
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
        case(3)
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
