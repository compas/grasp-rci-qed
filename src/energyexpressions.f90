!>     ./rci_qed_energyexpressions <state>
!!
!! Performs analysis and perturbative calculation on the RCI state.
!!
!! Takes a single command line argument: `state`.
!!
!! Depends on files: `<state>.c`
!!
program rci_qed_energyexpressions
    use grasp_rciqed_kinds, only: real64, dp
    use grasp_rciqed_lib9290_init, only: lib9290_init_constants, lib9290_init_rkco_gg
    use grasp_rciqed_lib9290_files, only: load_csl
    use grasp_rciqed_lib9290_csls, only: ncsfs_global

    use buffer_C, only: NVCOEF, COEFF, LABEL
    use orb_C, only: NW
    use prnt_C, only: NVEC
    use cord_I
    use onescalar_I
    use rkco_gg_I

    implicit none

    character(256) :: state
    integer :: state_len
    character(:), allocatable :: file_csls

    real(real64), allocatable :: tshell(:)
    integer :: i, j, k, l, m

    ! Fetch the first command line argument, which should be the name of the state.
    if (command_argument_count() /= 1) then
        error stop 'Missing command line argument: ./rci_qed_energyexpressions <STATE>'
    endif
    call get_command_argument(1, state) ! TODO: make state allocatable
    state_len = len_trim(state)

    ! The input files are assumed to be <state>.c, <state>.w and <state>.cm
    allocate(character(state_len+2) :: file_csls)
    file_csls = trim(state)//'.c'

    ! Make sure that all the files exist, or ERROR STOP otherwise
    call check_file(file_csls)

    ! Load the isotope information, CSL, radial orbitals and the CI mixing coefficents.
    print '(a)', ">>> INITIALIZING GLOBAL STATE"
    !call init_isocw_full(isodata, file_csls, file_wfns)
    call lib9290_init_constants ! physical and machine constants
    call load_csl(file_csls)
    call lib9290_init_rkco_gg ! TODO: does this need to be here?

    print '(a)', ">>> ENERGY EXPRESSIONS"
    allocate(tshell(NW))
    do i = 1, ncsfs_global()
        do j = 1, i
            call ONESCALAR(i, j, k, l, tshell)
            print '("ic, ir = ",i3,", ",i3)', i, j
            if(k /= 0 .and. k == l) then
                ! If k==l, then it seems that we have to add just the diagonal elements
                ! of the single particle matrix elements, weighted with TSHELL, to the
                ! matrix element.
                do m = 1, NW
                    print '(f10.5," (",i3,", ",i3,")")', tshell(m), m, m
                enddo
            elseif(k /= 0 .and. k /= l) then
                print '(f10.5," (",i3,", ",i3,")")', tshell(1), k, l
            endif

            ! Coulomb operator
            CALL RKCO_GG(i, j, CORD, 1, 1)
            do m = 1, NVCOEF
                print *, COEFF(m), LABEL(1,m), LABEL(2,m), LABEL(3,m), LABEL(4,m), LABEL(5,m)
            enddo
        enddo
    enddo

    print '(a)', ">>> Outputting data:"
    ! open(newunit=fh, file=trim(state)//".asfenergies.csv", action='write')
    ! call write_csv_header(fh)
    ! close(fh)

contains

    ! subroutine write_csv_header(unit)
    !     integer, intent(in) :: unit
    !
    !     write(unit, '(a9)', advance='no') "asf_index"
    !     do i = 1, hcs_len
    !         if(.not.hcs_cat(i)) cycle
    !         if(i == 7) then
    !             ! Self-energy gets special treatment
    !             write(unit, '(",",a24)', advance='no') trim(hcs_technical(i))//":"//trim(setypes_short(settings%qed_se))
    !         else
    !             write(unit, '(",",a24)', advance='no') trim(hcs_technical(i))
    !         endif
    !     enddo
    !     write(unit, '(",",a24)', advance='no') "sum:non_pt"
    !     do i = 1, hcs_len
    !         if(hcs_cat(i)) cycle
    !         if(i == 7) then
    !             ! Self-energy gets special treatment, since in PT mode, we have
    !             ! several parallel methods.
    !             do j = 1, nsetypes
    !                 write(unit, '(",",a24)', advance='no') "pt:"//trim(hcs_technical(i))//":"//trim(setypes_short(j))
    !             enddo
    !         else
    !             write(unit, '(",",a24)', advance='no') "pt:"//trim(hcs_technical(i))
    !         endif
    !     enddo
    !     write(unit, '()')
    ! end subroutine write_csv_header
    !
    ! subroutine write_csv_asf_line(unit, me)
    !     integer, intent(in) :: unit
    !     type(matrixelement), intent(in) :: me
    !
    !     real(real64) :: sum = 0.0_dp
    !
    !     write(unit, '(i9)', advance='no') k
    !     do i = 1, hcs_len
    !         if(.not.hcs_cat(i)) cycle
    !         write(unit, '(1(",",e24.16))', advance='no') getindex(me, i)
    !         sum = sum + getindex(me, i)
    !     enddo
    !     write(unit, '(1(",",e24.16))', advance='no') sum
    !     do i = 1, hcs_len
    !         if(hcs_cat(i)) cycle
    !         if(i == 7) then
    !             ! Self-energy gets special treatment, since in PT mode, we have
    !             ! several parallel methods.
    !             do j = 1, nsetypes
    !                 write(unit, '(1(",",e24.16))', advance='no') me%se_array(j)
    !             enddo
    !         else
    !             write(unit, '(1(",",e24.16))', advance='no') getindex(me, i)
    !         endif
    !     enddo
    !     write(unit, '()') ! write the newline
    ! end subroutine write_csv_asf_line

    subroutine check_file(filename)
        use grasp_rciqed_system, only: file_exists
        character(*), intent(in) :: filename

        if(.not.file_exists(filename)) then
            print '("Missing file ",a)', filename
            error stop
        endif
    end subroutine check_file

    ! !> Checks if the difference of `a` and `b` are within the tolerance relative
    ! !! to \f$\max(|a|,|b|)\f$.
    ! !!
    ! !! @param a,b Values to be checked.
    ! !! @param relative_tolerance Relative tolerance \f$\sigma\f$.
    ! !! @returns Whether \f$|a-b| / \max(|a|,|b|) < \sigma\f$.
    ! function within_tolerance(a, b, relative_tolerance)
    !     use grasp_rciqed_kinds, only: real64
    !
    !     real(real64), intent(in) :: a, b, relative_tolerance
    !     logical :: within_tolerance
    !     real(real64) :: relative_difference
    !
    !     relative_difference = abs(a-b) / max(abs(a), abs(b))
    !     if (relative_difference < relative_tolerance) then
    !         within_tolerance = .true.
    !     else
    !         within_tolerance = .false.
    !     endif
    ! end function within_tolerance
    !
    ! !> Calculate the relative difference of `a` and `b`.
    ! !!
    ! !! The relative difference is defined as:
    ! !!
    ! !! \f[
    ! !!   \frac{|a-b|}{\max(|a|, |b|)}
    ! !! \f]
    ! !!
    ! !! @param a,b Input values.
    ! !! @returns The relative difference of `a` and `b`.
    ! function reldiff(a, b)
    !     use grasp_rciqed_kinds, only: real64
    !     real(real64), intent(in) :: a, b
    !     real(real64) :: reldiff
    !     reldiff = abs(a-b) / max(abs(a), abs(b))
    ! end function reldiff
    !
    ! pure function check(yesno)
    !     logical, intent(in) :: yesno
    !     character(:), allocatable :: check
    !     if(yesno) then
    !         check = "✔"
    !     else
    !         check = "✘"
    !     endif
    ! end function check

end program rci_qed_energyexpressions
