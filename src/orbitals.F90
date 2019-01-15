!>     ./rci-qed.orbitals <state>
program rci_qed_orbitals
    use grasp_rciqed_system
    use grasp_rciqed_qed, only: init_vacuum_polarization
    use grasp_rciqed_qed_qedmod
    use grasp_lib9290, only: init_isocw_full
    use vacpol_I
    implicit none

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls, file_wfns

    character(*), parameter :: isodata = 'isodata' ! name of the isodata file

    ! Fetch the first command line argument, which should be the name of the state.
    if (command_argument_count() /= 1) then
        error stop 'Missing command line argument: ./rci_qedreport <STATE>'
    endif
    call get_command_argument(1, state)
    state_len = len_trim(state)

    ! The input files are assumed to be <state>.c and <state>.w
    file_csls = trim(state)//'.c'
    file_wfns = trim(state)//'.w'

    ! Make sure that all the files exist, or ERROR STOP otherwise
    call file_exists_or_stop(isodata,   __FILE__, __LINE__)
    call file_exists_or_stop(file_csls, __FILE__, __LINE__)
    call file_exists_or_stop(file_wfns, __FILE__, __LINE__)

    ! Load input files and initalize various parts of the GRASP global state
    print '(a)', ">>> Loading input files"
    call init_isocw_full(isodata, file_csls, file_wfns)

    ! Initialize the U+KS vacuum polarization potential
    call init_vacuum_polarization
    ! Initialize the qedmod common blocks with Grasp data etc.
    call qedse_qedmod_init

    call print_orbital_qed
    call qed_orbital_summary

contains

    !> Prints a table of QED self-energy estimates for all the orbitals.
    subroutine qed_orbital_summary
        use grasp_rciqed_kinds, only: real64, dp
        use grasp_rciqed_qed, only: qedse
        use orb_C, only: NW
        use vacpol_I
        implicit none

        real(real64) :: matrix(NW, NW)

        print *, 'QED operator matrix: vacuum polarization (U+KS)'
        call fill_vpint_matrix(matrix)
        call writematrix(matrix)

        print *, 'QED operator matrix: self-energy (Hydrogenic / Mohr)'
        call qedse(0, matrix)
        call writematrix(matrix)

        print *, 'QED operator matrix: self-energy (QEDMOD)'
        call qedse(1, matrix)
        call writematrix(matrix)

        print *, 'QED operator matrix: self-energy (Flambaum)'
        call qedse(2, matrix)
        call writematrix(matrix)

        print *, 'QED operator matrix: self-energy (Pyykkoe)'
        call qedse(3, matrix)
        call writematrix(matrix)
    end subroutine qed_orbital_summary

    subroutine print_orbital_qed
      use grasp_rciqed_kinds, only: real64, dp
      use grasp_rciqed_qed_pyykkoe
      use grasp_rciqed_qed_flambaum
      use grasp_rciqed_qed_qedmod
      use orb_C, only: NW, NP, NH
      use qed_slfen_I
      use vpint_I
      implicit none

      integer :: k
      real(real64) :: wfnorm, se_mohr(NW), se_pyykkoe, se_qedmod
      real(real64) :: se_flam, se_flam_phi_l, se_flam_phi_f, se_flam_phi_g
      real(real64) :: vp

      real(real64) :: se_qedmod_hydrogenic, qedse_qedmod_hydrogenic, dcwf_etot

      ! QED_SLFEN fills its argument array with Mohr & KLAMAQ self-energy values
      call QED_SLFEN(se_mohr)

      print *
      print *, 'QED self-energy values for orbitals'
      print '(a7,a12,9a15)', &
          '', 'WF', 'Mohr', &
          'Fl. low (l)', 'Fl. el (f)', 'Fl. mag (g)', 'Flambaum', &
          'Pyykkoe', 'QEDMOD', 'VP (U+KS)'

      do k = 1, NW
          wfnorm = qed_orbital_summary_wfnorm(k)
          se_flam = qedse_flambaum(k, k, se_flam_phi_l, se_flam_phi_f, se_flam_phi_g)
          se_pyykkoe = qedse_pyykkoe(k, k)
          se_qedmod = qedse_qedmod(k, k)
          call VPINT(k, k, vp)

          print '(i5,a2,es12.5,9es15.5)', &
              NP(k), NH(k), wfnorm, &
              se_mohr(k), &
              se_flam_phi_l, se_flam_phi_f, se_flam_phi_g, se_flam, &
              se_pyykkoe, se_qedmod, vp
      end do
    end subroutine print_orbital_qed

    subroutine fill_vpint_matrix(matrix)
        use grasp_rciqed_kinds, only: real64, dp
        use orb_C, only: NW, NAK
        use vpint_I
        implicit none

        real(real64), intent(out) :: matrix(NW, NW)

        real(real64) :: vp
        integer :: k, l

        do k = 1, NW
            do l = k, NW
                if(NAK(k) == NAK(l)) then
                    ! NOTE: VPINT swaps the indices if k > l, which will break the
                    ! do loop if that should happen. do l = k, NW ensures that this
                    ! does not happen.
                    call VPINT(k, l, vp)
                else
                    vp = 0_dp
                endif
                matrix(k, l) = vp
                matrix(l, k) = vp
            enddo
        enddo
    end subroutine fill_vpint_matrix

    subroutine writematrix(matrix)
        use grasp_rciqed_kinds, only: real64
        use orb_C, only: NW

        real(real64), intent(in) :: matrix(NW, NW)

        integer, allocatable :: kappas(:)
        integer :: idx

        kappas = get_kappas()
        do idx = 1, size(kappas)
            call writematrix_kappa(matrix, kappas(idx))
        enddo
    end subroutine writematrix

    subroutine writematrix_kappa(matrix, kappa)
        use grasp_rciqed_kinds, only: real64
        use orb_C, only: NW, NP, NH, NAK

        real(real64), intent(in) :: matrix(NW, NW)
        integer, intent(in) :: kappa

        integer :: k, l

        write(*, '(a5)', advance='no') ""
        do k = 1, NW
            if(NAK(k) /= kappa) then
                cycle
            endif
            write(*, '("           ",i2,a2)', advance='no') NP(k), NH(k)
        enddo
        write(*,'()')

        do k = 1, NW
            if(NAK(k) /= kappa) then
                cycle
            endif
            write(*, '(i2,a2)', advance='no') NP(k), NH(k)
            do l = 1, k-1
                if(NAK(l) /= kappa) then
                    cycle
                endif
                write(*, '(a15)', advance='no') ""
            enddo
            do l = k, NW
                if(NAK(l) /= kappa) then
                    cycle
                endif
                if(NAK(k) /= NAK(l)) then
                    ERROR STOP
                endif
                write(*, '(es15.5)', advance='no') matrix(k, l)
            enddo
            write(*,'()')
        enddo
    end subroutine writematrix_kappa

    !> Return an array of all the unique kappa values that the orbitals have.
    function get_kappas()
        use orb_C, only: NAK, NW

        integer, allocatable :: get_kappas(:)

        integer :: kappa, maxkappa, kappas_found = 0
        integer, allocatable :: tempkappas(:)

        maxkappa = MAXVAL(ABS(NAK(1:NW)))
        !print *, maxkappa ! TODO: clean this
        allocate(tempkappas(2*maxkappa))

        kappas_found = 0
        do kappa = 1, maxkappa
            if(has_kappa(-kappa)) then
                kappas_found = kappas_found + 1
                tempkappas(kappas_found) = -kappa
            endif

            if(has_kappa(kappa)) then
                kappas_found = kappas_found + 1
                tempkappas(kappas_found) = kappa
            endif
        enddo

        allocate(get_kappas(kappas_found))
        get_kappas(1:kappas_found) = tempkappas(1:kappas_found)
        deallocate(tempkappas)
    end function get_kappas
    !

    !> Is there a globally loaded GRASP orbital with the given `kappa` value?
    function has_kappa(kappa)
        use orb_C, only: NAK, NW

        integer, intent(in) :: kappa
        logical :: has_kappa

        integer :: k

        do k = 1, NW
            if(NAK(k) == kappa) then
                has_kappa = .true.
                return
            endif
        enddo
        has_kappa = .false.
    end function has_kappa

    function qed_orbital_summary_wfnorm(k)
        use grasp_rciqed_kinds, only: real64
        use grid_C, only: N, RP
        use wave_C, only: PF, QF
        use tatb_C, only: MTP, TA
        use quad_I
        implicit none

        integer :: k, i
        real(real64) :: qed_orbital_summary_wfnorm

        MTP = N
        do i = 1, N
            TA(i) = (PF(i,k)**2 + QF(i,k)**2) * RP(i)
        end do
        call QUAD(qed_orbital_summary_wfnorm)
    end function qed_orbital_summary_wfnorm

end program rci_qed_orbitals
