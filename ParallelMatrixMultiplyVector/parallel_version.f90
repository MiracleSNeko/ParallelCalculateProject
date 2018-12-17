program parallel_Mat_mul_Vec

    use mpi
    implicit none

    integer, parameter :: N = 1024
    integer :: IERR, NPROC, STATUS(MPI_STATUS_SIZE)
    integer :: myrank, myfile, buf_size
    real(4), allocatable :: matrix_buf(:, :), vector_buf(:) 

    call mpi_init(IERR)
    call mpi_comm_rank(MPI_COMM_WORLD, myrank, IERR)
    call mpi_comm_size(MPI_COMM_WORLD, NPROC, IERR)
    
    buf_size = N / NPROC
    ! read matrix
    call mpi_file_open(MPI_COMM_WORLD,"matrix",MPI_MODE_RDONLY &
    &  MPI_INFO_NULL, myfile, IERR)
    call mpi_file_seek(myfile, myrank*)
    call mpi_finalize(IERR)

end program parallel_Mat_mul_Vec