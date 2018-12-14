program generate_random_vector
    
    implicit none

    ! character :: s
    integer :: i
    integer, parameter :: n = 1024
    real(4), dimension(n) :: random_vector
    call random_seed()
    call random_number(random_vector)
    random_vector = 100 * random_vector

    open(20, file = 'vector.dat', status = 'replace')
    write(20, *) (random_vector(i), i = 1, n)
    close(20)

end program generate_random_vector

