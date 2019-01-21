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

    !> QED self-energy part of the `GETCID` routine.
    !!
    !! The questions get written to the standard error stream (`iounit_C :: istde`),
    !! in accordance with the rest of the `GETCID` implementation.
    subroutine getcid_qedse
        use iso_fortran_env, only: stdin => input_unit, stdout => output_unit
        use grasp_rciqed_cli, only: getoption
        use grasp_rciqed, only: setype
        use decide_C, only: LSE
        use qedcut_C, only: NQEDCUT, NQEDMAX
        use iounit_C, only: istde
        use getyn_I, only: getyn
        implicit none

        write(istde, *) "Estimate self-energy?"
        LSE = getyn()
        if(LSE .eqv. .true.) then
           NQEDCUT = 1
           print *, "Choose the method for QED self-energy estimation:"
           ! TODO: update the names
           write(istde, *) "    0 -- Hydrogenic / Mohr (default)"
           write(istde, *) "    1 -- QEDMOD (Shabaev et.al)"
           write(istde, *) "    2 -- Flambaum et.al"
           write(istde, *) "    3 -- Pyykkö et.al"
           setype = getoption(0, 3)
           if(setype == 0) then
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
