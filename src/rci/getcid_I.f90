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
        use grasp_rciqed, only: isspecorb

        integer :: specorbs(NNNW), specorbs_n

        write(istde, *) "List spectrosopic orbitals (e.g. for Breit, QED)"
        call GETRSL(specorbs, specorbs_n)
        call enabled_orbitals(specorbs(1:specorbs_n), isspecorb)
        ! TODO: print out a list of spectroscopic orbitals?
    end subroutine getcid_specorbs

    !> Converts a list of orbital indices `orbitals` into a list of `logical`s of length
    !! `NNNW`, where `.true.` values mark orbitals that are in `orbitals`.
    !!
    !! @param orbitals List of orbital indices.
    !! @param enabledorbitals Output array of `logical`s.
    !!
    !! Primarily a convenience function to convert the output of `GETRSL` into the input
    !! of `init_breit`.
    subroutine enabled_orbitals(orbitals, enabledorbitals)
        use parameter_def, only: NNNW
        integer, intent(in) :: orbitals(:)
        logical, intent(out) :: enabledorbitals(NNNW)
        integer :: i
        do i = 1, size(enabledorbitals)
            enabledorbitals(i) = .false.
        enddo
        do i = 1, size(orbitals)
            enabledorbitals(orbitals(i)) = .true.
        enddo
    end subroutine enabled_orbitals

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
        use grasp_rciqed, only: isspecorb, breit_mode
        use grasp_rciqed_breit, only: breit_specorbs

        character(len=1000) :: user_input
        integer :: i

        do while (.true.)
            write(istde, *) "Include contribution of H (Transverse)? Valid modes are:"
            write(istde, *) "   n        -- do not include H (Transverse)"
            write(istde, *) "   breit    -- enable in frequency-independent mode (Breit operator)"
            write(istde, *) "   specorbs -- enable full freq.dep. for spectroscopic orbitals, Breit for others"
            write(istde, *) "   full     -- enable full freq.dep. H (Transverse) for all orbitals (not recommended)"
            read(stdin, '(a)') user_input
            if(trim(user_input) == "n") then
                LTRANS = .false.
            elseif(trim(user_input) == "breit") then
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
