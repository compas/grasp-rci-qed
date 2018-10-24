program exportvp
    use ncdist_C, only: ZDIST
    use grid_C, only: N, R
    implicit none

    integer :: i, fh

    call setup
    call setup_vacuum_polarization

    open(newunit=fh, file="vacuum_polarization.csv", action='write')
    write(fh, '(a5,2(",",a25))') "idx", "radius", "vacuum_polarization"
    do i = 1, N
        write(fh, '(i5,2(",",es25.16))') i, R(i), ZDIST(i)
    enddo
    close(fh)

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

    subroutine setup_vacuum_polarization
        use decide_C, only: LVP
        use vpilst_C, only: FRSTVP, NVPI
        use ncdist_C, only: ZDIST
        use tatb_C, only: TB
        use grid_C, only: N, RP
        use ncharg_I
        use vacpol_I

        LVP = .true.
        call ncharg
        call vacpol
        ZDIST(2:N) = TB(2:N)*RP(2:N)
        FRSTVP = .TRUE.
        NVPI = 0
    end subroutine setup_vacuum_polarization

end program exportvp
