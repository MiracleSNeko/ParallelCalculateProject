!******************************************************************************
!
!  matrix_mul_vector_serial.f90
!  ���о��������˷�
!  
!******************************************************************************

program serial_Mat_mul_Vec

    implicit none

    integer, parameter :: N = 2048
    integer :: i, j
    real(4) :: startwtime, endwtime
    real(4) :: matrix(N, N), vector(N, 1), answer(N, 1) = 0
    
    call cpu_time(startwtime)

    ! ��ȡ����
    open(10, file = 'matrix', access = 'direct', form = 'unformatted', &
    &  recl = 4*N*N)
    read(10, rec = 1)  ((matrix(i,j), j = 1, N), i = 1, N)
    close(10)
    matrix = transpose(matrix)  ! Fortran�ľ���洢��ʽΪ�д洢����Ҫ
                                ! ����һ��ת��
    
    ! ��ȡ����
    open(20, file = 'vector', access = 'direct', form = 'unformatted', &
    &  recl = 4*N)
    read(20, rec = 1) (vector(i, 1), i = 1, N)
    close(20)
    
    answer = matmul(matrix, vector)
    
    ! �������������ļ�
    open(30, file = 'answer', access = 'direct', form = 'unformatted',&
    &  recl = 4)
    do i = 1, N
        write(30, rec = i) answer(i, 1)
    end do
    
    call cpu_time(endwtime)
    write(*, *) "walltime", ":", 1000 * (endwtime - startwtime), "ms"

end program serial_Mat_mul_Vec