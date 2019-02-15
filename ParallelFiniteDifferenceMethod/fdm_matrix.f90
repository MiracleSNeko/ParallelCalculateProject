!***********************************************************************
!
!        fdm_matrix.f90
!        ����ַ���Ⲵ�ɷ��̵�һ��ֵ����
!        ���̵ı�ֵ���������ݣ��Լ��������Ӱ��
!        \Delta u = xy  (x,y) \in (0,1) \times (0,1)
!        u(0,y) = y, u(1,y) = 1-y, u(x,0) = u(x,1) = 0
!        ��������־��󲢴����ֵ����
!
!***********************************************************************

!-------------------------������������----------------------------------

module constant_
    
    implicit none
    integer, parameter :: N = 50
    integer, parameter :: matrix_size = (N+1)*(N+1)
    real(4), parameter :: h = 1.0/N
    ! real(4) :: mat_A(matrix_size: matrix_size)
    ! real(4) :: vec_U(matrix_size: 1)
    ! real(4) :: vec_F(matrix_size: 1)
    integer :: i, j, k

end module constant_


!--------------------------���������-----------------------------------

program five_point_difference_matrix

    use constant_
    implicit none
    
    real(4) :: A(0: matrix_size-1, 0: matrix_size-1) = 0
    real(4) :: U(0: matrix_size-1, 1) = 0, F(0: matrix_size-1, 1) = 0
    

    

end program five_point_difference_matrix


!-------------------------�ӳ���ͺ�������------------------------------

subroutine generate_matrix(matrix_A, vector_F)

    use constant_
    implicit none
    
    real(4), intent(inout) :: matrix_A(0: matrix_size-1, 0: matrix_size-1)
    !real(4), intent(inout) :: vector_U(matrix_size, 1)
    real(4), intent(inout) :: vector_F(0: matrix_size-1, 1)
    
    ! ��������ֵ�ϵ������(�����߽�)
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


subroutine boundray_condition(matrix_A, vector_F)

    use constant_
    implicit none
    
    real(4), intent(inout) :: matrix_A(0: matrix_size-1, 0: matrix_size-1)
    real(4), intent(inout) :: vector_F(0: matrix_size-1, 1)
    
    ! ��ʼ���߽�����
    ! �ϱ߽�: U(x,1) = 0
    do j = 1, N
        matrix_A(N+j*(N+1), N+j*(N+1)) = 1
        vector_F(N+j*(N+1), 1) = 0
    end do
    ! �±߽�: U(x,0) = 0
    do j = 1, N
        matrix_A(0+j*(N+1), 0+j*(N+1)) = 1
        vector_F(0+j*(N+1), 1) = 0
    end do
    ! ��߽�: U(0,y) = y
    do i = 1,N
        matrix_A(i+0*(N+1), i+0*(N+1)) = 1
        vector_F(i+0*(N+1), 1) = h*i
    end do
    ! �ұ߽�: U(1,y) = 1-y
    do i = 1, N
        matrix_A(i+N*(N+1), i+N*(N+1)) = 1
        vector_F(i+N*(N+1), 1) = 1-h*i
    end do

    ! ���ϱ߽�����
    ! �����߽�
    
end subroutine boundray_condition

