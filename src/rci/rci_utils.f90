module grasp_rci_utils
    implicit none

contains

    subroutine write_zdist_csv
        use ncdist_C, only: ZDIST
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

end module grasp_rci_utils
