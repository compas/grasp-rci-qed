module matrix_I
    implicit none

    interface
        subroutine MATRIX(ncore, j2max)
            ! Implementation lives in matrix.f90
            integer, intent(in) :: ncore, j2max
        end subroutine MATRIX
    end interface

end module matrix_I
