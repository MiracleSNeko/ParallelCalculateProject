!***********************************************************************
!
!        fdm_matrix.f90
!        五点差分法求解泊松方程第一边值问题
!        方程的边值条件并不相容，对结果可能有一定影响
!        \Delta u = xy  (x,y) \in (0,1) \times (0,1)
!        u(0,y) = y, u(1,y) = 1-y, u(x,0) = u(x,1) = 0
!        生成五点差分矩阵并处理边值条件
!
!***********************************************************************

!-------------------------常数声明部分----------------------------------

module constant_
    
    implicit none
    integer, parameter :: N = 50
    integer, parameter :: matrix_size = (N+1)*(N+1)
    real(4), parameter :: h = 1.0/N
    integer, parameter :: num_not_bd = matrix_size-2*(N+1)-2*(N-1)
    integer :: i, j, row, col

end module constant_


!--------------------------主程序入口-----------------------------------

program five_point_difference_matrix

    use constant_
    implicit none
    
    real(4) :: A(0: matrix_size-1, 0: matrix_size-1) = 0
    real(4) :: F(0: matrix_size-1, 1) = 0
    integer :: color(num_not_bd)
    call generate_matrix(A, F)
    call boundray_condition(A, F)
    call sort_RB(A(0: num_not_bd-1, 0: num_not_bd-1), color)
    

end program five_point_difference_matrix


!-------------------------子程序与函数部分------------------------------



subroutine boundray_condition(matrix_A, vector_F)

    ! 处理完后仅保留A和F的一部分，其余部分置零

    use constant_
    implicit none
    
    real(4), intent(inout) :: matrix_A(0: matrix_size-1, 0: matrix_size-1)
    real(4), intent(inout) :: vector_F(0: matrix_size-1, 1)
    real(4), allocatable :: buf_A_1(:, :), buf_A_2(:, :), buf_F(:, :)
    integer :: index_zero(2*(N+1)-2*(N-1)), index_notzero(num_not_bd)
    
    ! 初始化边界条件
    ! 上边界: U(x,1) = 0
    j = N
    do i = 0, N+1
        matrix_A(i+j*(N+1), i+j*(N+1)) = 1
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 下边界: U(x,0) = 0
    j = 0
    do i = 0, N+1
        matrix_A(i+j*(N+1), i+j*(N+1)) = 1
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 左边界: U(0,y) = y
    i = 0
    do j = 0, N+1
        matrix_A(i+j*(N+1), i+j*(N+1)) = 1
        vector_F(i+j*(N+1), 1) = h*i
    end do
    ! 右边界: U(1,y) = 1-y
    i = N
    do j = 0, N+1
        matrix_A(i+j*(N+1), i+j*(N+1)) = 1
        vector_F(i+j*(N+1), 1) = 1-h*i
    end do

    ! 整合边界条件
    ! 各个边界减去对应列非零元所在行的元素
    ! 结果是将矩阵化为除边界外可以红黑排序的矩阵
    ! 处理上边界
    j = N
    do i = 0, N+1
        do col = 0, matrix_size
            if ( 0 /= matrix_A(i+j*(N+1), col))  then
                do row = 0, matrix_size
                    matrix_A(row, col) = matrix_A(row, col) - matrix_A(row, i+j*(N+1))
                end do
                vector_F(col, 1) = vector_F(col, 1) - vector_F(i+j*(N+1), 1)
            end if
        end do
        matrix_A(i+j*(N+1), i+j*(N+1)) = 0
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 处理下边界
    j = 0
    do i = 0, N+1
        do col = 0, matrix_size
            if ( 0 /= matrix_A(i+j*(N+1), col))  then
                do row = 0, matrix_size
                    matrix_A(row, col) = matrix_A(row, col) - matrix_A(row, i+j*(N+1))
                end do
                vector_F(col, 1) = vector_F(col, 1) - vector_F(i+j*(N+1), 1)
            end if
        end do
        matrix_A(i+j*(N+1), i+j*(N+1)) = 0
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 处理左边界
    i = 0
    do j = 0, N+1
        do col = 0, matrix_size
            if ( 0 /= matrix_A(i+j*(N+1), col))  then
                do row = 0, matrix_size
                    matrix_A(row, col) = matrix_A(row, col) - matrix_A(row, i+j*(N+1))
                end do
                vector_F(col, 1) = vector_F(col, 1) - vector_F(i+j*(N+1), 1)
            end if
        end do
        matrix_A(i+j*(N+1), i+j*(N+1)) = 0
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 处理右边界
    i = N
    do j = 0, N+1
        do col = 0, matrix_size
            if ( 0 /= matrix_A(i+j*(N+1), col))  then
                do row = 0, matrix_size
                    matrix_A(row, col) = matrix_A(row, col) - matrix_A(row, i+j*(N+1))
                end do
                vector_F(col, 1) = vector_F(col, 1) - vector_F(i+j*(N+1), 1)
            end if
        end do
        matrix_A(i+j*(N+1), i+j*(N+1)) = 0
        vector_F(i+j*(N+1), 1) = 0
    end do
    ! 去掉多余的零，得到可以红黑排序的矩阵
    ! 计算并记录全零行的位置
    col = 0
    ! 下边界
    j = 0
    do i = 0, N
        if (col >= 2*(N-1) + 2*(N+1)) exit
        index_zero(col) = i+j*(N+1)
        col = col + 1
    end do
    ! 上边界
    j = N
    do i = 0, N
        if (col >= 2*(N-1) + 2*(N+1)) exit
        index_zero(col) = i+j*(N+1)
        col = col + 1
    end do
    ! 左边界
    i = 0
    do j = 0, N
        if (col >= 2*(N-1) + 2*(N+1)) exit
        index_zero(col) = i+j*(N+1)
        col = col + 1
    end do
    ! 右边界
    i = N
    do j = 0, N
        if (col >= 2*(N-1) + 2*(N+1)) exit
        index_zero(col) = i+j*(N+1)
        col = col + 1
    end do
    ! 计算非全零行位置
    i = 0
    do j = 0, matrix_size
        do col = 0, 2*(N+1)+2*(N-1)
            if (j == index_zero(col)) exit
        end do
        if (col >= 2*(N+1)+2*(N-1)) then
            index_notzero(i) = j
            i = i + 1
        end if
    end do
    ! 去掉全零行
    allocate(buf_A_1(0:matrix_size-1, 0:num_not_bd-1))
    do i = 0, num_not_bd
        do j = 0 ,matrix_size
            buf_A_1(j, i) = matrix_A(j, index_notzero(i))
        end do
    end do
    ! 去掉全零列
    allocate(buf_A_2(0:num_not_bd-1, 0:num_not_bd-1))
    do i = 0, num_not_bd
        do j = 0 ,num_not_bd
            buf_A_1(j, i) = matrix_A(index_notzero(j), i)
        end do
    end do
    ! 对应处理右端项
    allocate(buf_F(0: num_not_bd-1, 1))
    do i = 0, num_not_bd
        buf_F(i, 1) = vector_F(index_notzero(i), 1)
    end do

    matrix_A = 0
    vector_F = 0
    matrix_A(0: num_not_bd-1, 0: num_not_bd-1) = buf_A_2(0: num_not_bd-1, 0: num_not_bd-1)
    vector_F(0: num_not_bd-1, 1) = buf_F(0: num_not_bd-1, 1)

    deallocate(buf_F)
    deallocate(buf_A_2)
    deallocate(buf_A_1)

end subroutine boundray_condition


subroutine sort_RB(matrix_A, c_)

    use constant_
    implicit none

    real(4), intent(in) :: matrix_A(num_not_bd, num_not_bd)
    integer :: l_ = 1 ! 颜色种数
    integer :: I_(num_not_bd) ! 所有节点组成的的集合
    integer, intent(inout) :: c_(num_not_bd) ! 各个节点所分配的颜色号
    integer :: T_(num_not_bd) ! 已标记过颜色的节点集合
    integer :: S_(num_not_bd) ! 尚未标记颜色的节点集合
    integer :: t(num_not_bd) ! 辅助数组
    integer :: adj_(num_not_bd) ! 邻节点组成的集合
    integer :: adj_cap_T(num_not_bd) ! 邻节点中已标记节点组成的集合
    integer :: k_ ! 节点k

    ! 各数组赋初值
    ! 对数组I_,T_,S_， 元素值为1表示对应下标的节点在数组中，0表示不在数组中
    ! 通过按位操作和Fortran隐式循环实现集合的运算
    c_(1) = l_
    I_ = 1
    T_ = 0
    S_ = 1
    c_(2: num_not_bd-1) = 0
    t_(2) = 0
    ! 遍历所有节点，对每个节点进行分类
    k_ = 0
    do while (0 == T_(num_not_bd))
        S_ = (/ (ieor(I_(i), T_(i)), i = 1, num_not_bd) /)
        ! 选取第一个尚未编号的节点
        do k_ = 1, num_not_bd
            if (1 == S_(k_)) exit
        end do
        ! 计算节点k_的全部邻节点
        do i = 1, num_not_bd
            if (0 /= matrix_A(k_, i)) adj_(i) = 1
        end do
        ! 计算邻节点中已标记部分
        adj_cap_T = (/ (iand(adj_(i), T_(i)), i = 1, num_not_bd) /)
        ! 更新辅助数组
        do j = 1, num_not_bd
            if (1 == adj_cap_T(j)) t(c_(j)) = 1
        end do
        ! 对未排序节点染色
        do j = 1, num_not_bd
            if (0 == t_(j)) exit
        end do
        c_(k_) = j
        ! 颜色数不足以完成排序，增加颜色
        if ( l_+1 == j) l_ = l_+1
        ! 更新数组T_
        do j = 1, num_not_bd
            if (1 == adj_cap_T(j)) t(c_(j)) = 0
        end do
        T_(k_) = ior(T_(k_), k_)
    end do

end subroutine sort_RB