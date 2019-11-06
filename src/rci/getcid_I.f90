module getcid_I
    implicit none

    interface
        subroutine GETCID(isofile, rwffile, idblk)
            ! Implementation lives in getcid.f90
            character(len=*) :: isofile, rwffile
            character(len=8), dimension(*) :: idblk
        end subroutine GETCID
    end interface

contains

    !> Part of the `GETCID` routine that reads in the list of spectrosopic orbitals from
    !! user.
    subroutine getcid_specorbs
        use parameter_def, ONLY: NNNW
        use getrsl_I
        use iounit_C, only: istde
        use grasp_rciqed, only: init_isspecorb

        integer :: specorbs(NNNW), specorbs_n

        write(istde, *) "List spectrosopic orbitals (e.g. for Breit, QED)"
        call GETRSL(specorbs, specorbs_n)
        print *, specorbs(1:specorbs_n) ! TODO: remove
        call init_isspecorb(specorbs_n, specorbs)
    end subroutine getcid_specorbs

    !> Breit interaction part of the `GETCID` routine.
    !!
    !! This populates two variables based on the user preferences:
    !!
    !! - `LTRANS` from `decide_C`: set to `.true.` if the user requested some form of Breit
    !!   to be included in the calculation, or `.false.` otherwise
    !! - `breit_specorbs` from `grasp_rciqed_breit`: if `LTRANS` is set, the orbitals for
    !!   which the full frequency-dependent Breit should be used, are set to `.true.` in
    !!   this array. For frequency-independent Breit, all elements are set to `.false.`.
    !!   If `LTRANS` is `.false.`, the contents of the array is undefined.
    !!
    subroutine getcid_breit
        use, intrinsic :: iso_fortran_env, only: stdin => input_unit
        use parameter_def, ONLY: NNNW
        use decide_C, only: LTRANS
        use iounit_C, only: istde
        use grasp_rciqed, only: isspecorb
        use grasp_rciqed_breit, only: breit_specorbs, breit_mode

        character(len=1000) :: user_input
        integer :: i

        do while (.true.)
            write(istde, *) "Include contribution of H (Transverse)? Valid modes are:"
            write(istde, *) "   n        -- do not include H (Transverse)"
            write(istde, *) "   gaunt    -- enable in frequency-independent mode"
            write(istde, *) "   specorbs -- enable full Breit only for spectroscopic orbitals, frequency-independent for others"
            write(istde, *) "   full     -- enable H (Transverse) for all orbitals (not recommended)"
            read(stdin, '(a)') user_input
            if(trim(user_input) == "n") then
                LTRANS = .false.
            elseif(trim(user_input) == "gaunt") then
                LTRANS = .true.
                breit_specorbs = (/(.false., i = 1, NNNW)/)
            elseif(trim(user_input) == "specorbs") then
                LTRANS = .true.
                breit_specorbs(:) = isspecorb(:)
            elseif(trim(user_input) == "full") then
                LTRANS = .true.
                breit_specorbs = (/(.true., i = 1, NNNW)/)
            else
                write(istde, *) "Invalid input for last option:", trim(user_input)
                cycle
            endif
            breit_mode = trim(user_input)
            exit
        enddo
    end subroutine getcid_breit

    !> QED self-energy part of the `GETCID` routine.
    !!
    !! The questions get written to the standard error stream (`iounit_C :: istde`),
    !! in accordance with the rest of the `GETCID` implementation.
    subroutine getcid_qedse
        use iso_fortran_env, only: stdin => input_unit
        use grasp_rciqed, only: setype
        use grasp_rciqed_cli, only: getoption
        use grasp_rciqed_qed, only: nsetypes, setypes_long
        use decide_C, only: LSE
        use qedcut_C, only: NQEDCUT, NQEDMAX
        use iounit_C, only: istde
        use getyn_I, only: getyn

        integer :: i

        write(istde, *) "Estimate self-energy?"
        LSE = getyn()
        if(LSE .eqv. .true.) then
            NQEDCUT = 1
            print *, "Choose the method for QED self-energy estimation:"
            do i = 1, nsetypes
                write(istde, '("   ",i2," -- ",a)') i, trim(setypes_long(i))
            enddo
            setype = getoption(1, nsetypes)
            if(setype == 1) then
                write(istde, *) "Largest n quantum number for including self-energy for orbital"
                write(istde, *) "n should be less or equal 8"
                read(stdin, *) NQEDMAX
            endif
        else
            setype = 0
            NQEDCUT = 0
        endif
    end subroutine getcid_qedse

end module getcid_I
