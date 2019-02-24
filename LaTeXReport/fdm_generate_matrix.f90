subroutine generate_matrix(matrix_A, vector_F)

    use constant_
    implicit none
    
    real(4), intent(inout) :: matrix_A(0: matrix_size-1, 0: matrix_size-1)
    real(4), intent(inout) :: vector_F(0: matrix_size-1, 1)
    
    `!` `构成五点差分的系数矩阵(不含边界)`
    do i = 1, N
        do j = 1, N
            matrix_A(i+j*(N+1), i+j*(N+1)) = -4
            matrix_A(i+(j-1)*(N+1), i+(j-1)*(N+1)) = 1
            matrix_A(i+(j+1)*(N+1), i+(j+1)*(N+1)) = 1
            matrix_A(i-1+j*(N+1), i-1+j*(N+1)) = 1
            matrix_A(i+1+j*(N+1), i+1+j*(N+1)) = 1
            vector_F(i+j*(N+1), 1) = i*j*h**4
        end do
    end do    

end subroutine generate_matrix