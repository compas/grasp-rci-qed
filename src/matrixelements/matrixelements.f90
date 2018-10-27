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
    use g2k_lib92
    use g2k_librci
    use g2k_csls
    implicit none

    interface
        function eval_matrixelement(ic, ir, hamcache)
            use g2k_librci
            integer, intent(in) :: ic, ir
            type(hamiltonian_cache), intent(in) :: hamcache
            type(matrixelement) :: eval_matrixelement
        end function eval_matrixelement
        subroutine eval_asfs(cimatrix)
            !use g2k_parameters
            use grasp_kinds, only: real64
            real(real64), dimension(:,:), intent(in) :: cimatrix
        end subroutine eval_asfs
    end interface

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls, file_wfns, file_mixing

    integer :: status, i, j
    integer :: fh

    type(hamiltonian_cache) :: hamcache
    type(matrixelement) :: hij

    real(real64), dimension(:,:), allocatable :: cimatrix

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
    print '(a)', ">>> CALLING: rcicommon_init"
    call rci_common_init

    ! Allocate the dense CI matrix
    allocate(cimatrix(ncsfs_global(), ncsfs_global()))

    ! Populate the Hamiltonian cache
    print '(a)', ">>> RUNNING: allocate_sematrices"
    call allocate_sematrices(hamcache)

    ! Diagonal matrix elements
    !call run_matrixelement(1, 1)

    open(newunit=fh, file="matrix.csv", action='write')
    write(fh, '(a6,",",a6)', advance='no') "row", "column"
    write(fh, '(4(",",a30))', advance='no') "diracpot", "coulomb", "breit", "vp"
    write(fh, '(4(",",a30))', advance='no') "se_mohr", "se_shabaev", "se_flambaum", "se_pyykkoe"
    write(fh, '(2(",",a30))', advance='no') "nms", "sms"
    write(fh, '()')
    do i = 1, ncsfs_global()
        do j = 1, i
            hij = eval_matrixelement(i, j, hamcache)
            call print_matrixelement(i, j, hij)
            call write_matrixelement(fh, i, j, hij)

            cimatrix(i, j) = hij%diracpot + hij%coulomb + hij%breit
            cimatrix(j, i) = cimatrix(i, j)
        enddo
    enddo
    close(fh)

    ! ...
    print '(a)', ">>> RUNNING: eval_asfs"
    call eval_asfs(cimatrix)

contains

    subroutine print_matrixelement(k, l, hij)
        use g2k_librci, only: matrixelement

        integer, intent(in) :: k, l
        type(matrixelement), intent(in) :: hij

        real(real64) :: total_dcb

        print '(80("="))'
        print '(">>> RUNNING: matrixelement(",i0,", ",i0,")")', k, l
        print '(80("-"))'
        print '("> Dirac + potential      = ",d20.10)', hij%diracpot
        print '("> Coulomb                = ",d20.10)', hij%coulomb
        print '("> Breit                  = ",d20.10)', hij%breit
        print '("> Vacuum polarization    = ",d20.10)', hij%vp
        print '("> Self-energy (Mohr)     = ",d20.10)', hij%se(1)
        print '("> Self-energy (Shabaev)  = ",d20.10)', hij%se(2)
        print '("> Self-energy (Flambaum) = ",d20.10)', hij%se(3)
        print '("> Self-energy (Pyykkö)   = ",d20.10)', hij%se(4)
        print '("> Normal mass shift      = ",d20.10)', hij%nms
        print '("> Special mass shift     = ",d20.10)', hij%sms

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
