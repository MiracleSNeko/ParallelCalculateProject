program generate_random_matrix
    
    implicit none

    ! character :: s
    integer :: i
    integer, parameter :: n = 1024
    real(4), dimension(n, n) :: random_matrix
    call random_seed()
    call random_number(random_matrix)
    random_matrix = 100 * random_matrix

    open(20, file = 'matrix', access = 'direct', &
    &  form = 'unformatted', recl = 4 * n)
    do i = 1, n
        write(20, rec = i) random_matrix(i, :)
    end do
    close(20)

end program generate_random_matrix

