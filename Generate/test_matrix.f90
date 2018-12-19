program generate_random_matrix
    
    implicit none

    ! character :: s
    integer :: i, j, k
    integer, parameter :: n = 32
    real(4) :: cnt
    
    open(20, file = 'testmatrix', access = 'direct', &
    &  form = 'unformatted', recl = 4)
    k = 0
    do i = 1, n
        do j = 1, n
            k = k + 1
            write(20, rec = k) cnt
            cnt = cnt + 1
        end do
    end do
    close(20)

end program generate_random_matrix

