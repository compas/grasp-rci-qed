!> Private module containing GRASP `lib92` common block definitions used by the
!! routines in `g2k_csls`.
!!
!! A separate module is used so that it would be possible to use `implicit none`
!! in the main function body.
!!
!! The module __should not__ be used outside of the `g2k_csls.f90` file.
!!
! module g2k_csls_common
!     IMPLICIT REAL*8 (A-H, O-Z)
!     include 'parameters.def'
!     COMMON/ORB2/NCF,NW,PNTRIQ &
!           /DEF2/C &
!           /DEF4/ACCY,NSCF,NSIC,NSOLV &
!           /DEF9/CVAC,PI &
!           /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
! end module g2k_csls_common

!> Routines used for interacting with the configuration state list (CSL) files.
module g2k_csls
    use orb_C
    use def_C
    use grid_C
    implicit none

contains

    !> Count the number of CSF blocks in the CSL file `filename`
    !!
    !! @param filename The name of the input configuration state list (CSL) file.
    !! @param nblocks The number of CSF blocks in `filename` (output).
    !! @param status Error status. `0` indicates a success, non-zero error (output).
    !!
    !! @todo The error status return value is currently not implemented.
    !! `count_blocks` simply `ERROR STOP`s when there is an error.
    subroutine count_blocks(filename, nblocks, status)
        character(*), intent(in) :: filename
        integer, intent(out) :: nblocks, status

        ! standard I/O error handling
        !  - iol stands for "I/O line (number)"
        ! The standard error handling label is 999
        integer :: ios, iol
        character(255) :: iom
        integer :: fhandle

        character(15) :: header_str
        character(2) :: linestart_str
        integer :: i

        ! TODO: Currently count_blocks only succeeds or ERROR STOPs.
        status = 0

        iol=__LINE__; open(newunit=fhandle, file=filename, form="formatted", status="old", iostat=ios, iomsg=iom, err=999)

        ! This is from the old SETCSLL routines
        iol=__LINE__; read(fhandle, '(1a15)', iostat=ios, iomsg=iom, err=999) header_str
        if(header_str /= "Core subshells:") then
            print *, "ERROR: Not a CSL file."
            print *, "Header does not match 'Core subshells:'"
            print *, "header_str:", header_str
            ERROR STOP ! TODO: Fix this -- decent error handling, please.
        endif

        ! Skip the next 4 records (from SETCSLL).
        ! Should skip the rest of the header.
        do i = 1, 4
            iol=__LINE__; read(fhandle, *)
            if(ios /= 0) then
                print *, "ERROR: Unable to read from", filename
                print *, ios, iom
                ERROR STOP ! TODO: Fix this -- decent error handling, please.
            endif
        enddo

        ! Start counting nblocks. The first CSL in the file does not start
        ! with a ' *', so we start counting from 1.
        nblocks = 1
        do
            ! if the following read fails, we assume that we have reached EOF
            read(fhandle,'(1a2)', iostat=ios, iomsg=iom) linestart_str
            if(ios /= 0) exit
            if(linestart_str == " *") then
                nblocks = nblocks + 1
            endif
        enddo

        close(fhandle)
        return

        ! I/O error handling
        999 continue
        print '(a)', "ERROR: IO error has occurred"
        print '(" in ",a,":",i0)', __FILE__, iol
        print '("iostat=", i0)', ios
        print '("iomsg=", a)', iom
        close(fhandle)
        ERROR STOP ! TODO: Fix this -- decent error handling, please.
    end subroutine count_blocks

    !> Returns the number of CSFs loaded to _global scope_.
    !!
    !! @return The number of CSFs loaded in global scope.
    function ncsfs_global()
        !use g2k_csls_common, only: NCF
        integer :: ncsfs_global
        ncsfs_global = NCF
    end function ncsfs_global

end module g2k_csls
