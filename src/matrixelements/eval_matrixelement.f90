subroutine eval_asfs(hamiltonian, k, hk)
    use grasp_kinds, only: real64, dp
    use g2k_librci, only: matrixelement
    use orb_C
    use prnt_C
    use syma_C
    use eigv_C
    use eigvec1_C
    implicit none

    type(matrixelement), dimension(:,:), intent(in) :: hamiltonian
    integer, intent(in) :: k
    type(matrixelement), intent(out) :: hk

    integer :: i, j
    real(real64) :: cicj

    print *, "NCF", NCF, "NW", NW
    !print *, NCFBL(1), NCFBL(2), NCFBL(3), NCFBL(4)
    print *, "NVEC", NVEC
    print *, "EAV", EAV
    print '(a10,2e24.15)', "EVAL(k)", EVAL(k), EVAL(k) + EAV
    do i = 1, NCF
        print *, EVEC(i + (k - 1) * NCF)
    enddo
    print *, "EVEC", EVEC(1), EVEC(2), EVEC(3)

    hk%diracpot = 0.0_dp
    hk%coulomb = 0.0_dp
    hk%breit = 0.0_dp
    hk%vp = 0.0_dp
    hk%nms = 0.0_dp
    hk%sms = 0.0_dp
    do i = 1, NCF
        do j = 1, NCF
            cicj = asf_coefficient(k, i) * asf_coefficient(k, j)
            hk%diracpot = hk%diracpot + hamiltonian(i, j)%diracpot * cicj
            hk%coulomb = hk%coulomb + hamiltonian(i, j)%coulomb * cicj
            hk%breit = hk%breit + hamiltonian(i, j)%breit * cicj
            hk%vp  = hk%vp  + hamiltonian(i, j)%vp * cicj
            hk%nms = hk%nms + hamiltonian(i, j)%nms * cicj
            hk%sms = hk%sms + hamiltonian(i, j)%sms * cicj
        enddo
    enddo

contains

    !> Returns the i-th ASF coefficient of the k-th state
    function asf_coefficient(k, i)
        integer, intent(in) :: k, i
        real(real64) :: asf_coefficient
        asf_coefficient = EVEC(i + (k - 1) * NCF)
    end function asf_coefficient

end subroutine eval_asfs
