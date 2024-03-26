module interface_cpp
    ! 声明外部C函数的接口
    interface
        subroutine serial_host(n, nz_num,ia, ja, a, x, w) bind(c)
            use iso_c_binding
            implicit none
            ! 接口声明中的参数类型和属性应与C函数的定义匹配
            integer(c_int), value :: n,nz_num
            ! type(c_ptr),value     :: ia,ja,a,x,w
            integer(c_int), dimension(n+1), intent(in) :: ia
            integer(c_int), dimension(nz_num), intent(in) :: ja
            real(c_double), dimension(nz_num), intent(in) :: a
            real(c_double), dimension(n), intent(in) :: x
            real(c_double), dimension(n), intent(inout) :: w
        end subroutine 
    end interface

    ! 在这里定义模块中的其他内容

end module interface_cpp