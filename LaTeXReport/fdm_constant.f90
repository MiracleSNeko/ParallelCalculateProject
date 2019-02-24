`!-------------------------常数声明部分----------------------------------`

module constant_
    
    implicit none
    integer, parameter :: N = 50
    integer, parameter :: matrix_size = (N+1)*(N+1)
    real(4), parameter :: h = 1.0/N
    integer, parameter :: num_not_bd = matrix_size-2*(N+1)-2*(N-1)
    integer :: i, j, row, col

end module constant_