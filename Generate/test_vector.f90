program generate_random_vector
    
    implicit none

    ! character :: s
    integer :: i
    integer, parameter :: n = 32
    real(4) :: cnt = 0
    

    open(20, file = 'testvector', access = 'direct', &
    &  form = 'unformatted', recl = 4)
    do i = 1, n
        write(20, rec = i) cnt
        cnt = 1 + cnt
    end do
    close(20)

end program generate_random_vector

