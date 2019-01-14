program generate_random_vector
    
    implicit none

    ! character :: s
    integer :: i
    integer, parameter :: n = 2048
    real(4), dimension(n) :: random_vector
    call random_seed()
    call random_number(random_vector)
    random_vector = 100 * random_vector

    open(20, file = 'vector', access = 'direct', &
    &  form = 'unformatted', recl = 4)
    do i = 1, n
        write(20, rec = i) random_vector(i)
    end do
    close(20)

end program generate_random_vector

