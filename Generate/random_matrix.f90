program generate_random_matrix
    
    implicit none

    ! character :: s
    integer :: i, j, k
    integer, parameter :: n = 1024
    real(4), dimension(n, n) :: random_matrix
    call random_seed()
    call random_number(random_matrix)
    random_matrix = 100 * random_matrix

    open(20, file = 'matrix', access = 'direct', &
    &  form = 'unformatted', recl = 4)
    k = 0
    do i = 1, n
        do j = 1, n
            k = k + 1
            write(20, rec = k) random_matrix(i, j)
        end do
    end do
    close(20)

end program generate_random_matrix

