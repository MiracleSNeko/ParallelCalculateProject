`!***********************************************************************`
`!`
`!`        `fdm\_matrix.f90`
`!`        `五点差分法求解泊松方程第一边值问题`
`!`        `方程的边值条件并不相容，对结果可能有一定影响`
`!`        `$\Delta u = xy,  (x,y) \in (0,1) \times (0,1)$`
`!`        `$u(0,y) = y, u(1,y) = 1-y, u(x,0) = u(x,1) = 0$`
`!`        `生成五点差分矩阵并处理边值条件`
`!`
`!***********************************************************************`


`!--------------------------主程序入口-----------------------------------`

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
