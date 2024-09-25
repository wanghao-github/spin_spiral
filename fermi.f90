module fermi_module
    implicit none
contains
    function fermi_array(omega, Beta_fake) result(value)
        ! 函数声明
        implicit none
        real(kind=8), intent(in) :: omega(:)        ! 输入数组
        real(kind=8), intent(in) :: Beta_fake       ! 标量参数
        real(kind=8), dimension(size(omega)) :: value   ! 返回数组
        integer :: i

        ! 计算费米分布函数
        do i = 1, size(omega)
            if (Beta_fake * omega(i) .ge. 20.0d0) then
                value(i) = 0.0d0
            elseif (Beta_fake * omega(i) .le. -20.0d0) then
                value(i) = 1.0d0
            else
                value(i) = 1.0d0 / (1.0d0 + exp(Beta_fake * omega(i)))
            endif
        end do
    end function fermi_array

    function fermi(omega, Beta_fake) result(value)
        implicit none
        real(kind(1.0d0)), intent(in) :: omega
        real(kind(1.0d0)), intent(in) :: Beta_fake
        real(kind(1.0d0)) :: value 
        if (Beta_fake*omega .ge. 20d0) then
            value = 0.0
        elseif (Beta_fake*omega .le. -20d0)then
            value = 1.0
        else
            value= 1.0/(1.0+exp(Beta_fake*omega))
   endif
   return
end function fermi    

    
end module fermi_module
