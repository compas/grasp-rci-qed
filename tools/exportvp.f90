program exportvp
    use, intrinsic :: iso_fortran_env, only: real64, dp => real64
    use grasp_rciqed_qed_vp, only: qedvp_init
    implicit none
    call setup
    call qedvp_init
    call write_zdist_csv

contains

    subroutine setup
        use parameter_def
        use def_C
        use grid_C
        use npar_C
        use tatb_C
        use nucpot_I
        use setcon_I
        use setmc_I
        use setqic_I
        use radgrd_I

        ! Set up constants
        ! ----------------------------------------------------------------------
        call setcon
        call setmc ! machine-dependent parameters


        ! Set up grid-related global state
        ! ----------------------------------------------------------------------
        H = 5.0D-2
        RNT = 2.0D-6
        HP = 0.0D0
        N = NNNP
        MTP = NNNP

        call setqic
        call radgrd

        ! Stolen from rci, based on the various examples of setting ACCY there
        ACCY = H**6


        ! Set up nucleus-related global state
        ! ----------------------------------------------------------------------

        ! C and CVAC are both speeds of light. However, C is usually read in from
        ! a file, so needs to be set manually.
        C = CVAC

        print *, TENMAX,EXPMAX,EXPMIN,PRECIS
        print *, CVAC, PI

        Z = 18.0D0
        NPARM = 2
        PARM(1) = 6.8839456274865651D-005
        PARM(2) = 9.8905913700962641D-006

        call nucpot
    end subroutine setup

    subroutine write_zdist_csv
        use grasp_rciqed_qed_vp, only: ZDIST
        use grid_C, only: N, R, RP
        implicit none

        integer :: i, fh

        open(newunit=fh, file="zdist.csv", action='write')
        write(fh, '(a5,3(",",a25))') "idx", "r", "rp", "zdist"
        do i = 1, N
            write(fh, '(i5,3(",",es25.16))') i, R(i), RP(i), ZDIST(i)
        enddo
        close(fh)
    end subroutine write_zdist_csv

end program exportvp
