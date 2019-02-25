subroutine sort_RB(matrix_A, c_)

    use constant_
    implicit none

    real(4), intent(in) :: matrix_A(num_not_bd, num_not_bd)
    integer :: l_ = 1 `!` `颜色种数`
    integer :: I_(num_not_bd) `!` `所有节点组成的的集合`
    integer, intent(inout) :: c_(num_not_bd) `!` `各个节点所分配的颜色号`
    integer :: T_(num_not_bd) `!` `已标记过颜色的节点集合`
    integer :: S_(num_not_bd) `!` `尚未标记颜色的节点集合`
    integer :: t(num_not_bd) `!` `辅助数组`
    integer :: adj_(num_not_bd) `!` `邻节点组成的集合`
    integer :: adj_cap_T(num_not_bd) `!` `邻节点中已标记节点组成的集合`
    integer :: k_ `!` `节点k`

    `!` `各数组赋初值`
    `!` `对数组I\_,T\_,S\_， 元素值为1表示对应下标的节点在数组中，0表示不在数组中`
    `!` `通过按位操作和Fortran隐式循环实现集合的运算`
    c_(1) = l_
    I_ = 1
    T_ = 0
    S_ = 1
    c_(2: num_not_bd-1) = 0
    t_(2) = 0
    `!` `遍历所有节点，对每个节点进行分类`
    k_ = 0
    do while (0 == T_(num_not_bd))
        S_ = (/ (ieor(I_(i), T_(i)), i = 1, num_not_bd) /)
        `!` `选取第一个尚未编号的节点`
        do k_ = 1, num_not_bd
            if (1 == S_(k_)) exit
        end do
        `!` `计算节点k\_的全部邻节点`
        do i = 1, num_not_bd
            if (0 /= matrix_A(k_, i)) adj_(i) = 1
        end do
        `!` `计算邻节点中已标记部分`
        adj_cap_T = (/ (iand(adj_(i), T_(i)), i = 1, num_not_bd) /)
        `!` `更新辅助数组`
        do j = 1, num_not_bd
            if (1 == adj_cap_T(j)) t(c_(j)) = 1
        end do
        `!` `对未排序节点染色`
        do j = 1, num_not_bd
            if (0 == t_(j)) exit
        end do
        c_(k_) = j
        `!` `颜色数不足以完成排序，增加颜色`
        if ( l_+1 == j) l_ = l_+1
        `!` `更新数组T\_`
        do j = 1, num_not_bd
            if (1 == adj_cap_T(j)) t(c_(j)) = 0
        end do
        T_(k_) = ior(T_(k_), k_)
    end do

end subroutine sort_RB