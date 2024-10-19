module orbital_op
    implicit none
contains

! 先产生一个num_wann*num_wann尺寸的0矩阵
! 判断左右 i j 的 m 来判断对应的矩阵元是否是零

! 写出矩阵元非零的规则

subroutine gen_orbital_matrix(num_wann,nrpts,hopping,orbital_momentum)
    ! complex,allocatable:: orbital_momentum(:,:,:,:)
    integer, intent(in) :: num_wann
    integer, intent(in) :: nrpts
    integer             :: ir, m, n
    complex             :: zi
    complex(kind(1.0d0)), intent(in)  :: hopping(num_wann,num_wann,nrpts)
    complex(kind(1.0d0)), intent(out) :: orbital_momentum(num_wann,num_wann,nrpts,3)
    orbital_momentum = (0,0)
    zi = (0.0,1.0)
    do ir=1, nrpts
        do m=1, num_wann
            if (orbital_nlm(m, 5) /= 1) cycle  ! 仅处理 l=1 的情况（p 轨道）
            do n = 1, num_wann
                if (orbital_nlm(n, 5) /= 1) cycle
                select case (orbital_nlm(m, 6))
                    case (-1)
                        if (orbital_nlm(n, 6) == 0) then
                            orbital_momentum(m, n, ir, 1) = -zi
                        elseif (orbital_nlm(n, 6) == 1) then
                            orbital_momentum(m, n, ir, 3) = zi
                        endif
                    case (0)
                        if (orbital_nlm(n, 6) == -1) then
                            orbital_momentum(m, n, ir, 1) = zi
                        elseif (orbital_nlm(n, 6) == 1) then
                            orbital_momentum(m, n, ir, 2) = -zi
                        endif
                    case (1)
                        if (orbital_nlm(n, 6) == 0) then
                            orbital_momentum(m, n, ir, 2) = zi
                        elseif (orbital_nlm(n, 6) == -1) then
                            orbital_momentum(m, n, ir, 3) = -zi
                        endif
                end select
            end do
        end do
    end do

    open(111,file='orbital_momentum')
        do dir =1,3
            do ir = 1, nrpts
                do m = 1, num_wann
                    do n= 1, num_wann
                        write(111,'(5I5,3F16.8)') irvec(1,ir) ,irvec(2,ir), irvec(3,ir), n,m, orbital_momentum(n,m,ir,dir)
                    enddo
                enddo
            enddo
        enddo
    close(111)