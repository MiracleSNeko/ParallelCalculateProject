!******************************************************************************
!
!  matrix_mul_vector_serial_with_mpi.f90
!  串行矩阵向量乘法
!  
!******************************************************************************

program serial_Mat_mul_Vec_mpiio

    use mpi
    implicit none

    integer, parameter :: N = 2048
    integer :: i, j
    integer :: IERR, NPROC, NSTATUS(MPI_STATUS_SIZE)
    integer :: myrank, myleft, myright, myfile, buf_size, cnt
    real(4) :: startwtime, endwtime, wtime
    real(4) :: matrix(N, N), vector(N, 1), answer(N, 1) = 0    
    
    call cpu_time(startwtime)
    
    call mpi_init(IERR)
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, IERR)
    call mpi_comm_size(MPI_COMM_WORLD, NPROC, IERR)
                                
    ! 读取矩阵
    call mpi_file_open(MPI_COMM_WORLD, "matrix", MPI_MODE_RDONLY, MPI_INFO_NULL, &
    &  myfile, IERR)
    call mpi_file_seek(myfile, myrank*N*N*sizeof(MPI_REAL), MPI_SEEK_SET, &
    &  IERR)
    call mpi_file_read(myfile, matrix, N*N, MPI_REAL, NSTATUS, IERR)
    call mpi_file_close(myfile, IERR)
    matrix = transpose(matrix)
    
    ! 读取向量
    call mpi_file_open(MPI_COMM_WORLD, "vector", MPI_MODE_RDONLY, MPI_INFO_NULL, &
    &  myfile, IERR)
    call mpi_file_seek(myfile, myrank*N*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_read(myfile, vector, N, MPI_REAL, NSTATUS, IERR)
    call mpi_file_close(myfile, IERR)

    
    answer = matmul(matrix, vector)
    
    ! 全规约结果向量，并行输出到文�?
    !call mpi_allreduce(answer, answer_buf, N, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, IERR)
    call mpi_file_open(MPI_COMM_WORLD, "answer", MPI_MODE_CREATE+MPI_MODE_WRONLY, &
    &  MPI_INFO_NULL, myfile, IERR)
    call mpi_file_seek(myfile, myrank*N*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_write(myfile, answer, N, MPI_REAL, &
    &  MPI_STATUS_IGNORE, IERR)
    call mpi_file_close(myfile, IERR)
    
    ! 将各进程的运行时间记录到文件�?
    call cpu_time(endwtime)
    wtime = (endwtime - startwtime) * 1000 
    call mpi_file_open(MPI_COMM_WORLD, "walltime_mpiio", MPI_MODE_CREATE &
    &  +MPI_MODE_WRONLY, MPI_INFO_NULL, myfile, IERR)
    call mpi_file_seek(myfile, myrank*sizeof(MPI_REAL), MPI_SEEK_SET, IERR)
    call mpi_file_write(myfile, wtime, 1, MPI_REAL, MPI_STATUS_IGNORE, IERR)
    call mpi_file_close(myfile, IERR)
    
    call mpi_finalize(IERR)

end program serial_Mat_mul_Vec_mpiio