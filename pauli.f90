module pauli_comp
    implicit none
    ! 定义双精度浮点数类型
    integer, parameter :: dp = selected_real_kind(15, 307)

contains

    subroutine pauli_block_all(M, size, pauli_result)

    !!! VASP wannier format up up dn dn
        implicit none
        integer, intent(in) :: size
        complex(kind(1.0d0)), intent(in) :: M(size, size)
        complex(kind(1.0d0)), intent(out) :: pauli_result(4)
        integer :: nwann_2
        complex(kind(1.0d0)) :: MI(size/2, size/2)
        complex(kind(1.0d0)) :: Mx(size/2, size/2)
        complex(kind(1.0d0)) :: My(size/2, size/2)
        complex(kind(1.0d0)) :: Mz(size/2, size/2)
        complex(kind(1.0d0)) :: trace_value

        nwann_2 = size / 2

        MI = (M(1:nwann_2, 1:nwann_2) + M(nwann_2+1:2*nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp
        Mx = (M(nwann_2+1:2*nwann_2, 1:nwann_2) + M(1:nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp
        My = (M(nwann_2+1:2*nwann_2, 1:nwann_2) - M(1:nwann_2, nwann_2+1:2*nwann_2)) * 0.5_dp
        Mz = (M(1:nwann_2, 1:nwann_2) - M(nwann_2+1:2*nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp

        call trace(MI, trace_value)
        pauli_result(1) = trace_value * 2.0_dp

        call trace(Mx, trace_value)
        pauli_result(2) = trace_value * 2.0_dp

        call trace(My, trace_value)
        pauli_result(3) = trace_value * 2.0_dp

        call trace(Mz, trace_value)

        pauli_result(4) = trace_value * 2.0_dp
    
    end subroutine pauli_block_all

subroutine pauli_block_all2(M, size, pauli_result)


        !!! QE wannier format up dn up dn
        implicit none
        integer, intent(in) :: size
        complex(kind(1.0d0)), intent(in) :: M(size, size)
        complex(kind(1.0d0)), intent(out) :: pauli_result(4)
        integer :: nwann_2
        complex(kind(1.0d0)) :: MI(size/2, size/2)
        complex(kind(1.0d0)) :: Mx(size/2, size/2)
        complex(kind(1.0d0)) :: My(size/2, size/2)
        complex(kind(1.0d0)) :: Mz(size/2, size/2)
        complex(kind(1.0d0)) :: trace_value

        nwann_2 = size / 2

        MI = (M(1:size:2, 1:size:2) + M(2:size:2, 2:size:2)) / 2.0_dp
        Mx = (M(1:size:2, 2:size:2) + M(2:size:2, 1:size:2)) / 2.0_dp
        My = (M(1:size:2, 2:size:2) - M(2:size:2, 1:size:2)) * 0.5_dp
        Mz = (M(1:size:2, 1:size:2) - M(2:size:2, 2:size:2)) / 2.0_dp

        call trace(MI, trace_value)
        pauli_result(1) = trace_value * 2.0_dp

        call trace(Mx, trace_value)
        pauli_result(2) = trace_value * 2.0_dp

        call trace(My, trace_value)
        pauli_result(3) = trace_value * 2.0_dp

        call trace(Mz, trace_value)
        pauli_result(4) = trace_value * 2.0_dp
    end subroutine pauli_block_all2


! def pauli_block_all(M):
!     MI = (M[::2, ::2] + M[1::2, 1::2]) / 2.0
!     Mx = (M[::2, 1::2] + M[1::2, ::2]) / 2.0
!     # Note that this is not element wise product with sigma_y but dot product
!     My = (M[::2, 1::2] - M[1::2, ::2]) * 0.5j
!     Mz = (M[::2, ::2] - M[1::2, 1::2]) / 2.0
!     return MI, Mx, My, Mz

    subroutine trace(M, trace_value)
        implicit none
        complex(kind(1.0d0)), intent(in) :: M(:,:)
        complex(kind(1.0d0)), intent(out) :: trace_value
        integer :: i, n

        n = size(M, 1)
        trace_value = 0.0_dp

        do i = 1, n
            trace_value = trace_value + M(i, i)
        end do
    end subroutine trace

end module pauli_comp
