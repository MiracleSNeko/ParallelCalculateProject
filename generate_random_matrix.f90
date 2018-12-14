program generate_random_matrix
    
    implicit none

    ! character :: s
    integer :: i
    integer, parameter :: n = 1024
    real(4), dimension(n, n) :: random_matrix
    call random_seed()
    call random_number(random_matrix)
    random_matrix = 100 * random_matrix

    open(20, file = 'matrix.dat', status = 'replace')
    write(20, *) (random_matrix(i, :), i = 1, n)
    close(20)

end program generate_random_matrix

