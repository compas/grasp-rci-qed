!>     ./rci.matrixelements <state>
!!
!! Performs analysis and perturbative calculation on the RCI state.
!!
!! Takes a single command line argument: `state`.
!!
!! Depends on files: `isodata`, `<state>.c`, `<state>.w` and `<state>.cm`
!!
program matrixelements
    use grasp_kinds, only: real64
    use grasp_rciqed_breit, only: init_breit
    use grasp_rciqed_mass_shifts, only: init_mass_shifts
    use grasp_rciqed_qed, only: init_vacuum_polarization
    use grasp_cimatrixelements
    use g2k_lib92
    use g2k_librci
    use g2k_csls
    use prnt_C, only: NVEC
    implicit none

    interface
        function eval_matrixelement(ic, ir, hamcache)
            use g2k_librci
            integer, intent(in) :: ic, ir
            type(hamiltonian_cache), intent(in) :: hamcache
            type(matrixelement) :: eval_matrixelement
        end function eval_matrixelement
        subroutine eval_asfs(hamiltonian, k, hk)
            use g2k_librci
            type(matrixelement), dimension(:,:), intent(in) :: hamiltonian
            integer, intent(in) :: k
            type(matrixelement), intent(out) :: hk
        end subroutine eval_asfs
    end interface

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls, file_wfns, file_mixing

    integer :: status, k, i, j, j2max
    integer :: fh

    type(hamiltonian_cache) :: hamcache
    type(matrixelement) :: hij

    type(matrixelement), dimension(:,:), allocatable :: hamiltonian
    real(real64) :: result
    real(real64), dimension(20) :: tshell

    character(*), parameter :: isodata = 'isodata' ! name of the isodata file

    ! Fetch the first command line argument, which should be the name of the state.
    if (command_argument_count() /= 1) then
        error stop 'Missing command line argument: ./rci_qedreport <STATE>'
    endif
    call get_command_argument(1, state)
    state_len = len_trim(state)

    ! The input files are assumed to be <state>.c, <state>.w and <state>.cm
    allocate(character(state_len+2) :: file_csls)
    file_csls = trim(state)//'.c'

    allocate(character(state_len+2) :: file_wfns)
    file_wfns = trim(state)//'.w'

    ! TODO: should be +3? Or rather, allocate is not actually necessary?
    allocate(character(state_len+2) :: file_mixing)
    file_mixing = trim(state)//'.cm'

    ! Make sure that all the files exist, or ERROR STOP otherwise
    call check_file(isodata)
    call check_file(file_csls)
    call check_file(file_wfns)
    call check_file(file_mixing)

    ! Load the isotope information, CSL, radial orbitals and the CI mixing coefficents.
    print '(a)', ">>> CALLING: lib92_init_cw"
    call lib92_init_cw(isodata, file_csls, file_wfns)
    print '(a)', ">>> CALLING: lib92_init_mixing"
    call lib92_init_mixing(file_mixing)
    print '(a)', ">>> INITIALIZING rci commons"
    call init_rkintc(j2max)
    call init_breit(j2max)
    call init_vacuum_polarization
    call init_mass_shifts

    ! Allocate the dense CI matrix
    allocate(hamiltonian(ncsfs_global(), ncsfs_global()))

    print '(a)', ">>> RUNNING: calculating matrixelements"
    open(newunit=fh, file="matrix.csv", action='write')
    write(fh, '(a6,",",a6)', advance='no') "row", "column"
    write(fh, '(4(",",a30))', advance='no') "diracpot", "coulomb", "breit", "vp"
    write(fh, '(4(",",a30))', advance='no') "se_mohr", "se_shabaev", "se_flambaum", "se_pyykkoe"
    write(fh, '(2(",",a30))', advance='no') "nms", "sms"
    write(fh, '()')
    do i = 1, ncsfs_global()
        do j = 1, i
            hamiltonian(i,j)%diracpot = dirac_potential(i, j)
            hamiltonian(i,j)%coulomb = coulomb(i, j)
            hamiltonian(i,j)%breit = breit(i, j)
            hamiltonian(i,j)%vp = qed_vp(i, j)
            hamiltonian(i,j)%nms = nms(i, j)
            hamiltonian(i,j)%sms = sms(i, j)
            call print_matrixelement(i, j, hamiltonian(i,j))
            if (i /= j) then
                hamiltonian(j,i) = hamiltonian(i,j)
            endif
        enddo
    enddo
    close(fh)

    ! ...
    print '(a)', ">>> RUNNING: eval_asfs"
    do k = 1, NVEC
        call eval_asfs(hamiltonian, k, hij)
        print '(a1,a10,a24)',    " ", "Type", "<H_i> (Ha)"
        print '(a1,a10,e24.15)', " ", "DC", hij%diracpot + hij%coulomb
        print '(a1,a10,e24.15)', " ", "Breit", hij%breit
        print '(a1,a10,e24.15)', ">", "DC+B", hij%diracpot + hij%coulomb + hij%breit
        print '(a1,a10,e24.15)', " ", "QED VP", hij%vp
        print '(a1,a10,e24.15)', " ", "NMS+SMS", hij%nms + hij%sms
        print '(a1,a10,e24.15)', ">", "Full", hij%diracpot + hij%coulomb + hij%breit &
            + hij%vp + hij%nms + hij%sms
    enddo

contains

    subroutine print_matrixelement(k, l, hij)
        use g2k_librci, only: matrixelement

        integer, intent(in) :: k, l
        type(matrixelement), intent(in) :: hij

        real(real64) :: total_dcb

        print '(80("="))'
        print '(">>> RUNNING: matrixelement(",i0,", ",i0,")")', k, l
        print '(80("-"))'
        print '("> Dirac + potential      = ",d24.15)', hij%diracpot
        print '("> Coulomb                = ",d24.15)', hij%coulomb
        print '("> Breit                  = ",d24.15)', hij%breit
        print '("> Vacuum polarization    = ",d24.15)', hij%vp
        print '("> Self-energy (Mohr)     = ",d24.15)', hij%se(1)
        print '("> Self-energy (Shabaev)  = ",d24.15)', hij%se(2)
        print '("> Self-energy (Flambaum) = ",d24.15)', hij%se(3)
        print '("> Self-energy (Pyykkö)   = ",d24.15)', hij%se(4)
        print '("> Normal mass shift      = ",d24.15)', hij%nms
        print '("> Special mass shift     = ",d24.15)', hij%sms

        total_dcb = hij%diracpot + hij%coulomb + hij%breit
        print '("H(",i0,", ",i0,"; DCB) = ",d20.10)', k, l, total_dcb
    end subroutine print_matrixelement

    subroutine write_matrixelement(fh, k, l, hij)
        use g2k_librci, only: matrixelement_t => matrixelement

        integer, intent(in) :: fh, k, l
        type(matrixelement), intent(in) :: hij

        write(fh, '(i6,",",i6)', advance='no') k, l
        write(fh, '(3(",",e30.20))', advance='no') hij%diracpot, hij%coulomb, hij%breit
        write(fh, '(1(",",e30.20))', advance='no') hij%vp
        write(fh, '(4(",",e30.20))', advance='no') hij%se(1), hij%se(2), hij%se(3), hij%se(4)
        write(fh, '(2(",",e30.20))', advance='no') hij%nms, hij%sms
        write(fh, '()') ! write the newline
    end subroutine write_matrixelement

    subroutine check_file(filename)
        character(*), intent(in) :: filename

        if(.not.file_exists(filename)) then
            print '("Missing file ",a)', filename
            error stop
        endif
    end subroutine check_file

    function file_exists(filename)
        character(*), intent(in) :: filename
        logical :: file_exists
        inquire(file=filename, exist=file_exists)
    end function file_exists

end program matrixelements
