!******************************************************************************
!
!  matrix_mul_vector_parallel.f90
!  ���о��������˷���������Ϊ�л��ַ������в��л�
!  
!******************************************************************************

program parallel_Mat_mul_Vec

    use mpi
    implicit none

    integer, parameter :: N = 2048
    integer :: my_left, my_right
    integer :: IERR, NPROC, NSTATUS(MPI_STATUS_SIZE)
    integer :: myrank, myleft, myright, myfile, buf_size, cnt
    real(4) :: startwtime, endwtime, wtime
    real(4), allocatable :: matrix(:, :), vector(:, :), answer(:, :) 
    real(4), allocatable :: matrix_buf(:, :), vector_buf(:, :), answer_buf(:, :)
    character(len = 2) :: sTemp
	
    call cpu_time(startwtime)
	
    call mpi_init(IERR)
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, IERR)
    call mpi_comm_size(MPI_COMM_WORLD, NPROC, IERR)
	   
    buf_size = N / NPROC
    myleft = my_left(myrank, NPROC)
    myright = my_right(myrank, NPROC)
	
    allocate(matrix_buf(N, buf_size)) ! Fortran�ľ��󴢴淽ʽΪ�д��棬��Ҫ
                                      ! ����һ��ת�ã�������ö�ȡ����ռ�
    allocate(vector_buf(buf_size,1)) ! �����������̴������������Ҫ�Ļ���ռ�
    allocate(matrix(buf_size, N))
    allocate(vector(buf_size, 1))
    allocate(answer(N, 1))
    allocate(answer_buf(N, 1))

    ! ��ȡ����
    call mpi_file_open(MPI_COMM_WORLD, "matrix", MPI_MODE_RDONLY, MPI_INFO_NULL, &
    &  myfile, IERR)
    call mpi_file_seek(myfile, myrank*N*buf_size*sizeof(MPI_REAL), MPI_SEEK_SET, &
    &  IERR)
    call mpi_file_read(myfile, matrix_buf, N*buf_size, MPI_REAL, NSTATUS, IERR)
    call mpi_file_close(myfile, IERR)
    matrix = transpose(matrix_buf)

    ! ��ȡ����
    call mpi_file_open(MPI_COMM_WORLD, "vector", MPI_MODE_RDONLY, MPI_INFO_NULL, &
    &  myfile, IERR)
    call mpi_file_seek(myfile, myrank*buf_size*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_read(myfile, vector, buf_size, MPI_REAL, NSTATUS, IERR)
    call mpi_file_close(myfile, IERR)

    answer = 0
    deallocate(matrix_buf) ! �ͷž��󻺴�ռ����ڴ���ÿһ�μ���ʱ�ľ����
    allocate(matrix_buf(buf_size, buf_size))
    ! ѭ�������д�������о����
    do cnt = 0, NPROC
        ! �����Ӧ������������ĳ˻�
        matrix_buf = matrix(:, mod(myrank+cnt, NPROC)*buf_size+1:(mod(myrank+cnt, NPROC) &
        &  +1)*buf_size)
        answer(myrank*buf_size+1:(myrank+1)*buf_size, :) = matmul(matrix_buf,vector) &
        &  + answer(myrank*buf_size+1:(myrank+1)*buf_size, :)
        ! ����һ��������Ĵ���(����)
        call mpi_send(vector, buf_size, MPI_REAL, myleft, myrank, MPI_COMM_WORLD, IERR)
        call mpi_recv(vector_buf, buf_size, MPI_REAL, myright, myright, &
        &  MPI_COMM_WORLD, NSTATUS, IERR)
        vector = vector_buf
    end do
	
    ! ȫ��Լ�������������������ļ�
    call mpi_allreduce(answer, answer_buf, N, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, IERR)
    call mpi_file_open(MPI_COMM_WORLD, "answer", MPI_MODE_CREATE+MPI_MODE_WRONLY, &
    &  MPI_INFO_NULL, myfile, IERR)
    call mpi_file_seek(myfile, myrank*buf_size*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_write(myfile, answer_buf(myrank*buf_size+1, 1), buf_size, MPI_REAL, &
    &  MPI_STATUS_IGNORE, IERR)
    call mpi_file_close(myfile, IERR)
	
    ! �������̵�����ʱ���¼���ļ���
    call cpu_time(endwtime)
    wtime = (endwtime - startwtime) * 1000
    write(sTemp, '(i2)') NPROC  
    call mpi_file_open(MPI_COMM_WORLD, "walltime"//trim(adjustl(sTemp)), MPI_MODE_CREATE &
    &  +MPI_MODE_WRONLY, MPI_INFO_NULL, myfile, IERR)
    call mpi_file_seek(myfile, myrank*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_write(myfile, wtime, 1, MPI_REAL, MPI_STATUS_IGNORE, IERR)
    call mpi_file_close(myfile, IERR)
	
    deallocate(matrix)
    deallocate(vector)
    deallocate(answer)
    deallocate(matrix_buf)
    deallocate(vector_buf)
    deallocate(answer_buf)
    call mpi_finalize(IERR)

end program parallel_Mat_mul_Vec


!-------------------�ӳ����뺯������-------------------------------------------

integer function my_left(myrank, nproc) result(ans)

    implicit none 
    integer, intent(in) :: myrank, nproc
    
    ans = myrank - 1
    if (0 == myrank) ans = nproc - 1

end function my_left


integer function my_right(myrank, nproc) result(ans)

    implicit none 
    integer, intent(in) :: myrank, nproc
    
    ans = myrank + 1
    if (nproc-1 == myrank) ans = 0

end function my_right